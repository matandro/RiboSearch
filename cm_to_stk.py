#!/usr/bin/env python3
'''
enter description
'''
import os
from typing import Dict, List, Tuple
from infernal import align_sequences


# Code to read sequences from FINAL_all.txt
def read_sequences(input_file: str) -> Dict[str, Dict[str, str]]:
    result = {}
    with open(input_file, 'r') as input_reader:
        input_reader.readline()
        for line in input_reader:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            seq_dict = result.get(parts[0], {})
            seq_dict[parts[1]] = parts[4]
            result[parts[0]] = seq_dict
    return result


# generates stockholm alignment files in sto_folder to all cm in cm_folder. assumes sequences in input_all
# Assumes cm files FINAL_<cm code>.cm where cm code is the identifier in input_all
def cm_to_sto(input_all: str, cm_folder: str, sto_folder: str):
    all_sequnces = read_sequences(input_all)
    cm_files = [file_name for file_name in os.listdir(cm_folder) if file_name[-3:] == '.cm'
                and file_name[:6] == 'FINAL_']
    for found_cm in cm_files:
        cm_code = found_cm[6:-3]
        tagged_sequences = all_sequnces.get(cm_code)
        if tagged_sequences is None:
            print("Could not file sequences for", found_cm)
            continue
        input_cm = os.path.join(cm_folder, found_cm)
        output_sto = os.path.join(sto_folder, 'FINAL_{}.sto'.format(cm_code))
        if not align_sequences(tagged_sequences, input_cm, output_sto, timeout=500):
            print("infernal.align_sequences failed to generate sto")


if __name__ == '__main__':
    BASE_FOLDER = "/opt/home/matan/Dropbox/PHd/RiboSearch/SandD_finals"
    cm_to_sto(os.path.join(BASE_FOLDER, "FINAL_all.txt"), BASE_FOLDER, BASE_FOLDER)
    BASE_FOLDER = "/opt/home/matan/Dropbox/PHd/RiboSearch/Dive Search/evalue0.01"
    cm_to_sto(os.path.join(BASE_FOLDER, "FINAL_all.txt"), BASE_FOLDER, BASE_FOLDER)
