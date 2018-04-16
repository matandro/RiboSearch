import logging
import sys
import os
from subprocess import Popen
from tempfile import NamedTemporaryFile as NTF
import vienna


VARNA_RNA = "/opt/algorithm/VARNArna/VARNAv3-93.jar"


def create_ct(sequence, structure):
    comp_map = {}
    closing_stack = []
    for i in range(0, len(sequence)):
        if structure[i] == '(':
            closing_stack.append(i)
        elif structure[i] == ')':
            start_index = closing_stack.pop()
            comp_map[start_index] = i
            comp_map[i] = start_index
    temp_file = NTF(dir = '.', delete = False, suffix='.ct', mode='w')
    for i in range(0, len(sequence)):
       next_item = i + 2
       if next_item > len(sequence):
           next_item = 0
       comp = comp_map.get(i)
       if comp is None:
           comp = 0
       else:
           comp += 1
       line = "{}\t{}\t{}\t{}\t{}\t{}\n".format(i + 1, sequence[i], i, next_item, comp, i) 
       temp_file.write(line)
    temp_file.close()
    return temp_file.name


def call_varna(sequence, structure, out_file_path):
    res = False
    ct_file = None 
    try:
        ct_file = create_ct(sequence, structure)
        param_list = ['java', '-cp', VARNA_RNA, 'fr.orsay.lri.varna.applications.VARNAcmd',
                  '-i', ct_file, '-o', out_file_path, '-resolution', '5.0']
        logging.info("Running {}".format(param_list))
        with Popen(param_list) as proc:
            proc.wait()
        res = True
    finally:
        if ct_file is not None and os.path.exists(ct_file):
            os.remove(ct_file)
    return res


def run_all_in_file(input_file_path, output_dir_path, seq_type='M'):
    if not os.path.isdir(output_dir_path):
        os.makedirs(output_dir_path)
    with open(input_file_path, 'r') as input_file:
        seq_index = 1
        for line in input_file:
            sequence = line.strip()
            if sequence == '':
                continue
            if seq_type == 'C':
                structure = vienna.fold(sequence)['centroid']
            else:
                structure = vienna.fold(sequence)['MFE']
            logging.info("Sequence {} run {}".format(seq_index, call_varna(sequence, structure, 
                                                                           os.path.join(output_dir_path, 
                                                                                        "seq_{}.jpg".format(seq_index)))))
            seq_index += 1


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    if len(sys.argv) < 4:
        logging.error('Usage: python3 generate_imgs.py <sequence list file> <image output directory> <M/C>')
        sys.exit(-1)
    run_all_in_file(sys.argv[1], sys.argv[2], sys.argv[3])

