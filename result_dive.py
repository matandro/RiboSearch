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


def get_tax_id(organism_id: str) -> int:
    req = requests.post('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',
                        data={'id': organism_id, 'db': 'nuccore', 'retmode': 'xml', 'rettype': 'fasta'})
    root = etree.fromstring(req.text)
    child_stack = [root]
    while child_stack:
        current = child_stack.pop()
        if current.tag == 'TSeq_taxid':
            return int(current.text)
        else:
            for child in current:
                child_stack.append(child)
    return None


def check_ancestor(name: str, tax_id: int, rank: str=None) -> bool:
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


def generate_clusters(match_file_path: str, design_file_path: str, is_filter_bacteria: bool=False) \
        -> List[DesignGroup]:
    design_group_map = {}
    with open(match_file_path, 'r') as match_file, open(design_file_path, 'r') as design_file:
        seq_code_map = {}
        design_file.readline()
        for line in design_file:
            items = line.strip().split('\t')
            seq_code_map[items[0]] = {'sequence': items[2].strip(), 'structure': items[4].strip()}
        match_file.readline()
        for line in match_file:
            items = line.strip().split('\t')
            design_id = items[0].strip()
            design_group = design_group_map.get(design_id, DesignGroup(design_id, seq_code_map.get(design_id)))
            if not is_filter_bacteria or not check_ancestor('Bacteria', get_tax_id(items[5].strip().split('/', 1)[0])):
                design_group.add_match(items[5].strip(), {'identifier': items[5].strip(), 'sequence': items[2].strip(),
                                                          'round': 0})
            design_group_map[design_id] = design_group
    return list(design_group_map.values())


def dive_single(group_id: str, single_design_group: DesignGroup, cm_dir: str, seq_db_path: str,
                filter_evalue: float = 10.0) -> Tuple[DesignGroup, int]:
    count = 0
    found_new = True
    base_cm_name = '{}.cm'.format(group_id)
    cm_name = 'TEMP_{}'.format(base_cm_name)
    shutil.copyfile(os.path.join(cm_dir, base_cm_name), os.path.join(cm_dir, cm_name))
    stockholm_file = '{}.stk'.format(group_id)
    design_group_identifies = {'sequence': single_design_group.sequence, 'structure': single_design_group.structure}
    design_copy = copy(single_design_group)
    if not os.path.exists(os.path.join(cm_dir, cm_name)):
        full_list = {}
        for identifier, match in single_design_group.matches.items():
            full_list[identifier] = match.get('sequence')
        full_list[single_design_group.identifier] = single_design_group.sequence
        success = infernal.align_sequences(full_list,
                                           os.path.join(cm_dir, cm_name), stockholm_file)
        success = infernal.generate_cm(stockholm_file, os.path.join(cm_dir, cm_name))
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
        success = infernal.generate_cm(stockholm_file, os.path.join(cm_dir, cm_name))
        # search on cm
        search_res = infernal.search_cm(os.path.join(cm_dir, cm_name), seq_db_path)
        # identify items (see different matches) and compare size of match group
        new_design_group = DesignGroup(single_design_group.identifier, design_group_identifies)
        for single_match in search_res:
            code = single_match.get('identifier')
            seq = single_match.get('sequence')
            if float(single_match.get('E-value')) < filter_evalue:
                new_design_group.add_match(code, single_match)
                if code not in design_copy.matches:
                    single_match['round'] = count
                    found_new = True
        design_copy = new_design_group
    # organize cm
    shutil.move(os.path.join(cm_dir, cm_name), os.path.join(cm_dir, 'FINAL_{}'.format(base_cm_name)))
    return design_copy, count


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


def run_dive(base_dir: str, filter_evalue: float = 10.0, filter_path: str = None):
    def check_filter():
        if filter_path is not None:
            with open(filter_path, 'r') as filter_file:
                filter_list = [line.strip() for line in filter_file if line.strip() != '']
                if not filter_list:
                    logging.error("filter file empty {}".format(filter_path))
                    exit(-1)

    logging.basicConfig(level=logging.INFO)
    filter_list = None
    logging.info('Reading clusters {} and {}'.format(os.path.join(base_dir, 'match_log'),
                                                     os.path.join(base_dir, 'design_log')))
    all_design_groups = generate_clusters(os.path.join(base_dir, 'match_log'), os.path.join(base_dir, 'design_log'))
    logging.info('Read {} clusters, starting dive'.format(len(all_design_groups)))
    check_filter()
    with open(os.path.join(base_dir, 'FINAL_summary'), 'w') as out_file, \
            open(os.path.join(base_dir, 'FINAL_all.txt'), 'w') as final_all:
        final_all.write('design_code\tidentifier\tscore\tE-value\tsequence\tRound\n')
        out_file.write('design code\tOriginal # of matches\tDive # of matches\t# of cycles\thas non bacteria\tsequence'
                       '\tstructure\n')
        for design_group in all_design_groups:
            if filter_list is not None and \
                    design_group.identifier not in filter_list:
                continue
            logging.info('Dive group {}'.format(design_group.identifier))
            new_design_group, count = dive_single(design_group.identifier, design_group, base_dir, '/DB/fasta_db/nt/nt',
                                                  filter_evalue)
            logging.info('Finished group {} matched {} results'.format(design_group.identifier,
                                                                       len(new_design_group.matches)))
            tax_name_list = [x.split('/')[0] for x in new_design_group.matches.keys()]
            write_str = '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(new_design_group.identifier, len(design_group.matches),
                                                            len(new_design_group.matches), count,
                                                            has_non_bacteria(tax_name_list), new_design_group.sequence,
                                                            new_design_group.structure)
            for identifier, match in new_design_group.matches.items():
                final_all.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(new_design_group.identifier, identifier,
                                                              match.get('score'), match.get('E-value'),
                                                              match.get('sequence'), match.get('round')))
            logging.info(write_str)
            out_file.write('{}\n'.format(write_str))
    logging.info('All clusters done')


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    # test_taxonomy()
    if len(sys.argv) == 1:
        base_dir = '/DB/Output/SandD'
        logging.info("No arguments, running on default ({}, 0.01, no filter)".format(base_dir))
        run_dive(base_dir, 0.01)
    elif len(sys.argv) < 3:
        logging.error("If arguments are given 2 are mandatory. result_dive.py <base_dir> <filter_amount> [filter path]"
                      "\nGot: {}".format(sys.argv))
        sys.exit(-1)
    else:
        logging.info("running with arguments: {}".format(sys.argv))
        run_dive(sys.argv[1], float(sys.argv[2]), None if len(sys.argv) < 4 else sys.argv[3])
