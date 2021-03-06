#!/usr/bin/env python3
'''
For each results, rebuild model and keep searching. Use infernal cm align to build a new model.
'''

import requests
from typing import List, Dict, Tuple
from ete3 import NCBITaxa
import xml.etree.ElementTree as etree
import infernal
import os
import sys
import shutil
from copy import copy
import logging
from rnafbinv import vienna, shapiro_tree_aligner
import dive_statistics
import argparse
from enum import Enum
import time


# workaround to the issue of search_runner having two different versions
class MODE(Enum):
    OLD = 1
    NEW = 2


FILTER_THRESHOLD = 300
mapped_map = {}

def get_tax_id(organism_id: str, max_retry: int=5) -> int:
    result = mapped_map.get(organism_id)
    if result is not None:
        return result
    works = False
    attempt = 0
    while not works:
        attempt += 1
        try:
            req = requests.post('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',
                                data={'id': organism_id, 'db': 'nuccore', 'retmode': 'xml', 'rettype': 'fasta'})
            if req.status_code < 200 or req.status_code > 300:
                raise requests.exceptions.RequestException("Bad response code {}".format(req.status_code))
            root = etree.fromstring(req.text)
            works = True
            child_stack = [root]
            while child_stack:
                current = child_stack.pop()
                if current.tag == 'TSeq_taxid':
                    mapped_map[organism_id] = int(current.text)
                    return int(current.text)
                else:
                    for child in current:
                        child_stack.append(child)
        except requests.exceptions.RequestException as re:
            print("get_tax_id request exception for organism {}\n{}".format(organism_id, re))
            if attempt > max_retry:
                break
            else:
                time.sleep(5)
        except Exception as e:
            print("FETAL ERROR FOR {}\n{}".format(organism_id, e))
    return None


def check_ancestor(name: str, tax_id: int, rank: str = None) -> bool:
    ncbi = NCBITaxa()
    ancestor_ids = ncbi.get_name_translator([name]).get(name, [])
    if not ancestor_ids:
        raise ValueError("No taxonomy id for {}".format(name))
    lineage = ncbi.get_lineage(tax_id)
    for anc_id in lineage:
        if rank is None or ncbi.get_rank([anc_id]).get(anc_id, '') == rank:
            if anc_id in ancestor_ids:
                return True
    return False


class DesignGroup:
    def __init__(self, code: str, seq_code_map: Dict[str, str]):
        self.identifier = code
        self.sequence = seq_code_map.get('sequence')
        self.structure = seq_code_map.get('structure')
        self.matches = {}

    def add_match(self, identifier: str, infernal_datum: Dict[str, str]):
        datum = self.matches.get(identifier)
        if datum is not None and datum != infernal_datum:
            logging.warning('Identifier [{}] already exists with different sequence:\n{}\n{}'.format(identifier,
                                                                                                     infernal_datum,
                                                                                                     datum))
        self.matches[identifier] = infernal_datum


def generate_clusters(match_file_path: str, design_file_path: str, target_tree,
                      is_filter_bacteria: bool = False, mode: MODE = MODE.OLD) -> List[DesignGroup]:
    vienna_folder = None
    try:
        if mode == MODE.NEW:
            vienna_folder = vienna.LiveRNAfold()
            vienna_folder.start()
        design_group_map = {}
        with open(match_file_path, 'r') as match_file, open(design_file_path, 'r') as design_file:
            seq_code_map = {}
            design_file.readline()
            for line in design_file:
                if line.strip() == '':
                    continue
                items = line.strip().split('\t')
                if mode == MODE.OLD:
                    seq_code_map[items[0]] = {'sequence': items[2].strip(), 'structure': items[4].strip()}
                else:
                    # new doesnt go through a filter so we will filter it now
                    code = '{}_{}'.format(items[0].strip(), items[1].strip())
                    sequence = items[3].strip()
                    structure = vienna_folder.fold(sequence)['MFE']
                    source_tree = shapiro_tree_aligner.get_tree(structure, sequence)
                    _, score = shapiro_tree_aligner.align_trees(source_tree, target_tree)
                    if score < FILTER_THRESHOLD:
                        seq_code_map[code] = {'sequence': sequence, 'structure': structure}
            match_file.readline()
            for line in match_file:
                if line.strip() == '':
                    continue
                items = line.strip().split('\t')
                design_id = items[0].strip()
                if seq_code_map.get(design_id) is None:
                    continue
                design_group = design_group_map.get(design_id, DesignGroup(design_id, seq_code_map.get(design_id)))
                if mode == MODE.OLD:
                    if not is_filter_bacteria or not check_ancestor('Bacteria',
                                                                    get_tax_id(items[5].strip().split('/', 1)[0])):
                        design_group.add_match(items[5].strip(), {'identifier': items[5].strip(),
                                                                  'sequence': items[2].strip(), 'round': 0})
                else:
                    # new mode didnt save identifier, add just and replace on first search
                    design_group.add_match(str(len(design_group.matches)), {'identifier': len(design_group.matches),
                                                                            'sequence': items[1].strip(), 'round': 0})
                design_group_map[design_id] = design_group
    finally:
        if vienna_folder is not None:
            vienna_folder.close()
    return list(design_group_map.values())


def get_align_score(identifier: str, sequence: str, cm_path: str, target_tree):
    fasta_file = None
    try:
        fasta_file = infernal.generate_fasta({identifier: sequence})
        cm_struct, new_sequence = dive_statistics.get_cm_struct(cm_path, fasta_file.name)
        cm_tree = shapiro_tree_aligner.get_tree(cm_struct, new_sequence)
        _, score = shapiro_tree_aligner.align_trees(cm_tree, target_tree)
        return score
    finally:
        if fasta_file is not None:
            os.remove(fasta_file.name)


def dive_single(group_id: str, single_design_group: DesignGroup, cm_dir: str, seq_db_path: str, target_tree,
                filter_align_score: float = 250, filter_evalue: float = 10.0, cpus: int = 12) -> \
        Tuple[DesignGroup, int, Dict[int, List[str]]]:
    count = 0
    items_in_round = {}
    found_new = True
    base_cm_name = '{}.cm'.format(group_id)
    if not os.path.exists(os.path.join(cm_dir, base_cm_name)):
        infernal.generate_single_seq_cm(single_design_group.sequence, os.path.join(cm_dir, base_cm_name),
                                        structure=single_design_group.structure, cpus=cpus)
    cm_name = 'TEMP_{}'.format(base_cm_name)
    shutil.copyfile(os.path.join(cm_dir, base_cm_name), os.path.join(cm_dir, cm_name))
    stockholm_file = os.path.join(cm_dir, '{}.sto'.format(group_id))
    design_group_identifies = {'sequence': single_design_group.sequence, 'structure': single_design_group.structure}
    design_copy = copy(single_design_group)
    items_in_round[0] = design_copy.matches.keys()
    while found_new:
        count += 1
        found_new = False
        # rebuild cm (align to old, delete and create new)
        full_list = {}
        for identifier, match in single_design_group.matches.items():
            full_list[identifier] = match.get('sequence')
        full_list[single_design_group.identifier] = single_design_group.sequence
        success = infernal.align_sequences(full_list,
                                           os.path.join(cm_dir, cm_name), stockholm_file)
        os.remove(os.path.join(cm_dir, cm_name))
        cm_path = os.path.join(cm_dir, cm_name)
        success = infernal.generate_cm(stockholm_file, cm_path, cpus=cpus)
        # search on cm
        search_res = infernal.search_cm(cm_path, seq_db_path, cpus=cpus)
        # identify items (see different matches) and compare size of match group
        new_design_group = DesignGroup(single_design_group.identifier, design_group_identifies)
        for single_match in search_res:
            code = single_match.get('identifier')
            seq = single_match.get('sequence')
            align_score = get_align_score(code, seq, cm_path, target_tree)
            if float(single_match.get('E-value')) < filter_evalue and align_score < filter_align_score:
                old_res = design_copy.matches.get(code)
                if old_res is None:
                    single_match['round'] = count
                    found_new = True
                else:
                    single_match['round'] = old_res['round']
                new_design_group.add_match(code, single_match)
        design_copy = new_design_group
        items_in_round[count] = design_copy.matches.keys()
    # organize cm
    shutil.move(os.path.join(cm_dir, cm_name), os.path.join(cm_dir, 'FINAL_{}'.format(base_cm_name)))
    shutil.move(stockholm_file, os.path.join(cm_dir, 'FINAL_{}.sto'.format(group_id)))
    return design_copy, count, items_in_round


def test_taxonomy():
    # Test Bacteria
    pulled_tax_id = get_tax_id('AE017333.1')
    print("TAX_ID", pulled_tax_id)
    print("Is bacteria?", check_ancestor('Bacteria', pulled_tax_id))
    # Test non Bacteria
    pulled_tax_id = get_tax_id('XM_014337919.1')
    print("TAX_ID", pulled_tax_id)
    print("Is bacteria?", check_ancestor('Bacteria', pulled_tax_id))


def has_non_bacteria(tax_list: List[str]) -> bool:
    has_warning = False
    for tax_name in tax_list:
        try:
            if not check_ancestor('Bacteria', get_tax_id(tax_name)):
                return True
        except:
            has_warning = True
            logging.warning("Unknown taxonomy {}".format(tax_name))
    if has_warning:
        return None
    else:
        return False


def run_dive(base_dir: str, target_tree, filter_align_score: float = 250, filter_evalue: float = 10.0,
             filter_path: str = None, mode: MODE = MODE.OLD, cpus: int = None):
    def check_filter():
        show_list = []
        if filter_path is not None:
            with open(filter_path, 'r') as filter_file:
                show_list = [line.strip() for line in filter_file if line.strip() != '']
                if not show_list:
                    logging.error("filter file empty {}".format(filter_path))
                    exit(-1)
        return show_list

    logging.basicConfig(level=logging.INFO)
    logging.info("Running run_dive with: Dir: {}, Max evalue: {}, Filter file: {}".format(base_dir, filter_evalue,
                                                                                          filter_path))
    logging.info('Reading clusters {} and {}'.format(os.path.join(base_dir, 'match_log'),
                                                     os.path.join(base_dir, 'design_log')))
    all_design_groups = generate_clusters(os.path.join(base_dir, 'match_log'), os.path.join(base_dir, 'design_log'),
                                          target_tree, mode=mode)
    logging.info('Read {} clusters, starting dive'.format(len(all_design_groups)))
    filter_list = check_filter()
    existed = os.path.exists(os.path.join(base_dir, 'FINAL_summary'))
    with open(os.path.join(base_dir, 'FINAL_summary'), 'a+') as out_file, \
            open(os.path.join(base_dir, 'FINAL_all.txt'), 'a+') as final_all, \
            open(os.path.join(base_dir, 'FINAL_per_round'), 'a+') as final_round:
        if not existed:
            final_all.write('design_code\tidentifier\tscore\tE-value\tsequence\tRound\n')
            final_all.flush()
            out_file.write(
                'design_code\tOriginal # of matches\tDive # of matches\t# of cycles\thas non bacteria\tsequence'
                '\tstructure\n')
            out_file.flush()
            final_round.write('design_code\tround\t# of items\titem identifiers\n')
            final_round.flush()
        for design_group in all_design_groups:
            if filter_list is not None and \
                    design_group.identifier not in filter_list:
                logging.info('+++++ {} not in filter file, continuing'.format(design_group.identifier))
                continue
            final_cm_path = os.path.join(base_dir, 'FINAL_{}.cm'.format(design_group.identifier))
            if os.path.exists(final_cm_path) and os.stat(final_cm_path).st_size > 0:
                logging.info('+++++ Already done  group {}'.format(design_group.identifier))
                continue
            logging.info('+++++ Dive group {}'.format(design_group.identifier))
            new_design_group, count, items_per_round = dive_single(design_group.identifier, design_group,
                                                                   base_dir, '/DB/fasta_db/nt/nt', target_tree,
                                                                   filter_evalue=filter_evalue, cpus=cpus,
                                                                   filter_align_score=filter_align_score)
            logging.info('Finished group {} matched {} results'.format(design_group.identifier,
                                                                       len(new_design_group.matches)))
            tax_name_list = [x.split('/')[0] for x in new_design_group.matches.keys()]
            write_str = '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(new_design_group.identifier, len(design_group.matches),
                                                            len(new_design_group.matches), count,
                                                            has_non_bacteria(tax_name_list), new_design_group.sequence,
                                                            new_design_group.structure)
            for round, ident_list in items_per_round.items():
                final_round.write('{}\t{}\t{}\t{}\n'.format(new_design_group.identifier, round, len(ident_list),
                                                            ", ".join(ident_list)))
            for identifier, match in new_design_group.matches.items():
                final_all.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(new_design_group.identifier, identifier,
                                                                  match.get('score'), match.get('E-value'),
                                                                  match.get('sequence'), match.get('round')))
            logging.info(write_str)
            out_file.write('{}\n'.format(write_str))
            out_file.flush()
            final_round.flush()
            final_all.flush()
    logging.info('All clusters done')


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    # test_taxonomy()
    if len(sys.argv) == 1:
        target_tree = shapiro_tree_aligner.get_tree(
            "((((((((...(.(((((.......))))).)........((((((.......))))))..))))))))",
            "NNNNNNNNUNNNNNNNNNNNNNNNNNNNNNNNNUNNNUNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNN")
        base_dir = '/DB/Output/SandD'
        logging.info("No arguments, running on default ({}, 0.01, no filter)".format(base_dir))
        run_dive(base_dir, target_tree, filter_evalue=0.01, mode=MODE)
    else:
        parser = argparse.ArgumentParser(description="A tool that performs family generation to search_runner.py single"
                                                     " results", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("base_dir", help="Path to directory containing the match_log and design_log files")
        parser.add_argument("-e", "--evalue", type=float, help="cmsearch evalue cutoff for results", default=10.0)
        parser.add_argument("-m", "--mode", type=str, help="modes: new for new search runner, old otherwise",
                            default='old', choices=['old', 'new'])
        parser.add_argument("-s", "--score", type=float, help="score cuttoff based on RNAfbinv alignment score",
                            default=250)
        parser.add_argument("-t", "--tree", type=str, help="sequence and structure seperated by _ of design targets",
                            default="NNNNNNNNUNNNNNNNNNNNNNNNNNNNNNNNNUNNNUNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNN_"
                                    "((((((((...(.(((((.......))))).)........((((((.......))))))..))))))))")
        parser.add_argument("-n", type=int, help="Work on design target with at least n minimal matches", default=5)
        parser.add_argument("--cpu", type=int, help="maximum number of CPU to use for infernal programs")
        parser.add_argument("-f", "--filter", type=str, help="path to file with newline separated design id's to work "
                                                             "on. ignores all other designs")
        args = parser.parse_args()
        seq_struct = args.tree.strip('"').strip("'").split('_')
        target_tree = shapiro_tree_aligner.get_tree(seq_struct[1].strip(), seq_struct[0].strip())
        run_dive(args.base_dir, target_tree, filter_evalue=args.evalue, filter_align_score=args.score,
                 cpus=args.cpu, filter_path=args.filter, mode=MODE[args.mode.upper()])
