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


def get_tax_id(organism_id: str) -> str:
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
    def __init(self, seq_code_map: Dict[str, str]):
        self.sequence = seq_code_map.get('sequence')
        self.structure = seq_code_map.get('structure')
        self.matches = []

    def add_match(self, sequence: str, structure:str=None):
        self.matches.append(sequence)


def generate_clusters(match_file_path: str, design_file_path: str, is_filter_bacteria: bool=False) \
        -> Dict[str, DesignGroup]:
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
            design_group = design_group_map.get(design_id, DesignGroup(seq_code_map.get(design_id)))
            if not is_filter_bacteria or not check_ancestor('Bacteria', get_tax_id(items[5].strip().split('/', 1)[0])):
                design_group.add_match(items[2].strip(), items[3].strip())
    return design_group_map


def dive_single(group_id: str, single_design_group: DesignGroup, cm_dir: str, seq_db: str):
    found_new = True
    base_cm_name = '{}.cm'.format(group_id)
    cm_name = base_cm_name
    stockholm_file = '{}.stk'.format(group_id)
    round = 1
    while found_new:
        found_new = False
        # rebuild cm
        success = infernal.align_sequences(dict(enumerate(single_design_group.matches)), os.path.join(cm_dir, cm_name)
                                           , stockholm_file)
        cm_name = '{}_{}'.format(round, base_cm_name)
        success = infernal.generate_cm(stockholm_file, os.path.join(cm_dir, cm_name))
        # search with new cm
        search_res = infernal.search_cm(os.path.join(cm_dir, cm_name), seq_db)
        # identify items (see different matches) and compare size of match group
        for single_match in search_res:
            code = '{}:{}-{}({})'.format(single_match.get('target name'), single_match.get('seq from'),
                                         single_match.get('seq to'), single_match.get('strand'))

    # organize cm
'''
target name': res_name, 'score': res_info[3], 'E-value': res_info[2], 
                   'seq from': res_info[9], 'seq to': res_info[10], 'strand': res_info[11],
                   'model from': res_info[6], 'model to': res_info[7], 'sequence': res_seq}
    '''



if __name__ == "__main__":
    # Test Bacteria
    pulled_tax_id = get_tax_id('AE017333.1')
    print("TAX_ID", pulled_tax_id)
    print("Is bacteria?", check_ancestor('Bacteria', pulled_tax_id))
    # Test non Bacteria
    pulled_tax_id = get_tax_id('XM_014337919.1')
    print("TAX_ID", pulled_tax_id)
    print("Is bacteria?", check_ancestor('Bacteria', pulled_tax_id))