#!/usr/bin/env python3
'''
For each results, rebuild model and keep searching. Use infernal cm align to build a new model.
'''

import requests
from typing import List, Dict
from ete3 import NCBITaxa
import xml.etree.ElementTree as etree
import infernal
import os
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

    def add_match(self, identifier: str, sequence: str, structure:str=None):
        seq = self.matches.get(identifier)
        if seq is not None and seq != sequence:
            logging.warning('Identifier [{}] already exists with different sequence:\n{}\n{}'.format(identifier, seq,
                                                                                                     sequence))
        self.matches[identifier] = sequence


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
                design_group.add_match(items[5].strip(), items[2].strip(), items[3].strip())
            design_group_map[design_id] = design_group
    return list(design_group_map.values())


def dive_single(group_id: str, single_design_group: DesignGroup, cm_dir: str, seq_db_path: str) -> DesignGroup:
    found_new = True
    base_cm_name = '{}.cm'.format(group_id)
    cm_name = 'TEMP_{}'.format(base_cm_name)
    shutil.copyfile(os.path.join(cm_dir, base_cm_name), os.path.join(cm_dir, cm_name))
    stockholm_file = '{}.stk'.format(group_id)
    design_group_identifies = {'sequence': single_design_group.sequence, 'structure': single_design_group.structure}
    design_copy = copy(single_design_group)
    while found_new:
        found_new = False
        # rebuild cm (align to old, delete and create new)
        success = infernal.align_sequences(single_design_group.matches,
                                           os.path.join(cm_dir, cm_name), stockholm_file)
        os.remove(os.path.join(cm_dir, cm_name))
        success = infernal.generate_cm(stockholm_file, os.path.join(cm_dir, cm_name))
        # search on cm
        search_res = infernal.search_cm(os.path.join(cm_dir, cm_name), seq_db_path)
        # identify items (see different matches) and compare size of match group
        newl_design_group = DesignGroup(single_design_group.identifier, design_group_identifies)
        for single_match in search_res:
            code = single_match.get('identifier')
            seq = single_match.get('sequence')
            newl_design_group.add_match(code, seq)
            if code not in design_copy.matches:
                found_new = True
        design_copy = newl_design_group
    # organize cm
    shutil.move(os.path.join(cm_dir, cm_name), os.path.join(cm_dir, 'FINAL_{}'.format(base_cm_name)))
    return design_copy


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


if __name__ == "__main__":
    # test_taxonomy()
    logging.basicConfig(level=logging.INFO)
    base_dir = '/DB/Output/SandD'
    logging.info('Reading clusters {} and {}'.format(os.path.join(base_dir, 'match_log'),
                                                     os.path.join(base_dir, 'design_log')))
    all_design_groups = generate_clusters(os.path.join(base_dir, 'match_log'), os.path.join(base_dir, 'design_log'))
    logging.info('Read {} clusters, starting dive'.format(len(all_design_groups)))
    with open(os.path.join(base_dir, 'FINAL_summary'), 'w') as out_file:
        out_file.write('design code\tNo of matches\thas non bacteria\tsequence\tstructure\n')
        for design_group in all_design_groups:
            logging.info('Dive group {}'.format(design_group.identifier))
            new_design_group = dive_single(design_group.identifier, design_group, base_dir, '/DB/fasta_db/nt/nt')
            logging.info('Finished group {} matched {} results'.format(design_group.identifier,
                                                                       len(new_design_group.matches)))
            tax_name_list = [x.split('/')[0] for x in new_design_group.matches.keys()]
            write_str = '{}\t{}\t{}\t{}\t{}'.format(new_design_group.identifier, len(new_design_group.matches),
                                                    has_non_bacteria(tax_name_list),
                                                    new_design_group.sequence, new_design_group.structure)
            logging.info(write_str)
            out_file.write('{}\n'.format(write_str))
    logging.info('All clusters done')
