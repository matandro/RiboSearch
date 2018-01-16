from subprocess import Popen, PIPE
import logging
import sys


RNAINV_DISTANCE_PATH = "/home/matan/Dropbox/Thesis/Workspace/RNAinv/RNAfbinv/RNAshapiroSeq"


#0> Enter sequence or structure:
#uuggaa
#1> seq = uuggaa
#2> energy score = 0.000000
#3> Fold = ......
#4> b2C = (R)
#5> shapiro = ((E6)R)
#6> shapiro clean = (R)
#7> other seq = AAGGUU
#8> other energy score = 0.000000
#9> other fold = ......
#10> other b2C = (R)
#11> other shapiro = ((E6)R)
#12> other shapiro clean = (R)
#13> shapiro distance = 0
#14> bp distance = 0
def _read_output_single(output):
    logging.debug("Reading output: \n[" + output + "]")
    result = {}
    if output is not None and "" != output:
        try:
            lines = output.split('\n')
            fold_energy = float(lines[2].split('=')[1].strip())
            result['energy'] = fold_energy
            result['structure'] = lines[3].split('=')[1].strip()
            result['shapiro'] = lines[5].split('=')[1].strip()
        except:
            pass

    if len(result) == 0:
        result = None
        logging.debug("Could not read output")
    return result


def _read_output_cmp(output):
    result = _read_output_single(output)
    if result is not None:
        try:
            lines = output.split('\n')
            fold_energy = float(lines[8].split('=')[1].strip())
            result['energy_cmp'] = fold_energy
            result['structure_cmp'] = lines[9].split('=')[1].strip()
            result['shapiro_cmp'] = lines[11].split('=')[1].strip()
            result['shapiro_distance'] = int(lines[13].split('=')[1].strip())
            result['bp_distance'] = int(lines[14].split('=')[1].strip())
        except:
            pass
    return result


def get_fold(sequence):
    return _run_base([RNAINV_DISTANCE_PATH], [sequence], _read_output_single)


def get_fold_compare(sequence_one, sequence_cmp):
    return _run_base([RNAINV_DISTANCE_PATH, sequence_cmp], [sequence_one], _read_output_cmp)


def _run_base(param_list, write_list, output_read_func):
    try:
        with Popen(param_list, stdout=PIPE, stdin=PIPE, stderr=PIPE) as proc:
            logging.debug("Running get_fold on params: {}".format(param_list))
            logging.debug("writing: [{}]".format('\n'.join(write_list)))
            fold_output = proc.communicate(input='\n'.join(write_list).encode())[0]
            logging.debug('output received: [{}]'.format(fold_output))
            result = output_read_func(fold_output.decode())
    except:
        e = sys.exc_info()[0]
        logging.error("Failed to run: '{}'\nerror: {}".format(param_list, e))
        sys.exc_info()
        sys.exit(1)
    return result


if __name__ == "__main__":
    #DEBUG
    #logging.basicConfig(level= logging.DEBUG)
    #argv = ['rnainv_folder', 'GUCACGUACGACGCGGAUAAUACCAUACGUGUUUACGACGAUC']
    # REAL
    argv = sys.argv
    argc = len(argv)
    if argc == 2:
        print(get_fold(argv[1]))
    elif argc == 3:
        print(get_fold_compare(argv[1], argv[2]))
    else:
        logging.error("Usage: python rnainv_folder <fasta sequence> [fasta sequence cmp]")
