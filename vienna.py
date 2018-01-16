from subprocess import Popen, PIPE
import logging
import os
import sys


RNAFOLD_PATH = "/opt/algorithm/ViennaRNA/bin/"
RNAFOLD_EXE = "RNAfold"


def output_fold_analyze(output):
    structure = {}
    try :
        lines = output.split('\\n')
        structure['MFE'] = lines[1].split(' ')[0].strip()
        structure['centroid'] = lines[3].split(' ')[0].strip()
    except:
        structure = {}
    return structure


def fold(sequence):
    structure_map = None
    try:
        param_list = [os.path.join(RNAFOLD_PATH, RNAFOLD_EXE), '-p']
        with Popen(param_list, stdout=PIPE, stdin=PIPE) as proc:
            logging.debug("Running RNAfold: {}".format(param_list))
            fold_output, errs = proc.communicate("{}\n@\n".format(sequence).encode())
            structure_map = output_fold_analyze(str(fold_output))
    except OSError as e:
        logging.error("Failed to run: '{}'. ERROR: {}".format(param_list, e.errno))
    return structure_map


if __name__ == "__main__":
    test_sequence = "AGUAGAUGGCCCGUGGUGUCCCGGAGUGGCUGUAGAGUGAGAUGCAGAUGGAC" \
                    "GACUGAGCCCAUAGGGCCGCUUAUAAUAACCUAUGCCCCCACAUCGUGUAAUU" \
                    "UCAACCCGCAGCACUAUCACAGCCACAGGGUCGAUCA"
    if len(sys.argv) > 1:
        test_sequence = sys.argv[1]
    print(fold(test_sequence))
