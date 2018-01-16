from tempfile import NamedTemporaryFile as NTF
from subprocess import Popen, PIPE, CalledProcessError
import logging
import sys
import os
import random


#RNAFBINV_PATH = '/home/matan/Dropbox/Thesis/Workspace/RNAinv/RNAfbinv/RNAexinv'
RNAFBINV_PATH = '/opt/algorithm/RNAinv/RNAfbinv/RNAexinv'


def _handle_output(output):
    logging.debug("Reading output: \n[" + output + "]")
    result = None
    if output is not None and "" != output:
        try:
            lines = output.split('\n')
            for line in lines:
                if line.startswith('Output sequence'):
                    result = line.split('=')[1].strip()
        except:
            pass
    if result is None:
        logging.error("Could not read output [{}]".format(output))
    return result


def _run_base(param_list):
    try:
        with Popen(param_list, stdout=PIPE, stdin=PIPE, stderr=PIPE) as proc:
            logging.info("Running RNAexinv on params: {}".format(param_list))
            rnafbinv_output, rnafbinv_error = proc.communicate()
            error_text = rnafbinv_error.decode()
            if error_text.strip() != '':
                logging.error("Error from RNAexinv: {}".format(error_text))
                raise CalledProcessError()
            result = _handle_output(rnafbinv_output.decode())
            logging.debug('exit code: [{}]'.format(proc.wait()))
    except:
        e = sys.exc_info()
        logging.error("Failed to run: '{}'\nerror: {}".format(param_list, e))
        sys.exc_info()
        sys.exit(1)
    return result


def run_rnafbinv(structure, sequence, energy=-1000, mr=-1000, iter_no=500, seed=None, random_seed=None):
    tmp_file = NTF(dir='.', delete=False)
    tmp_file.write('{}\n'.format(structure).encode())
    tmp_file.write('{}\n'.format(sequence).encode())
    tmp_file.write('{}\n'.format(energy).encode())
    tmp_file.write('{}\n'.format(mr).encode())
    tmp_file.close()
    try:
        param_list = [RNAFBINV_PATH, '-f', tmp_file.name]
        if seed is not None:
            param_list.extend(['-S', seed])
        if iter_no is not None:
            param_list.extend(['-i', str(iter_no)])
        if random_seed is None:
            random_seed = random.getrandbits(64)
        param_list.extend(['-p', str(random_seed)])
        result = _run_base(param_list)
    finally:
        os.remove(tmp_file.name)
    return result


if __name__ == "__main__":
    # DEBUG
    #logging.basicConfig(level=logging.DEBUG)
    argv = sys.argv
    if len(sys.argv) < 3:
        argv = ['','((((((((...(.(((((.......))))).)........((((((.......))))))..))))))))',
                'NNNNNNNNUNNNNNNNNNNNNNNNNNNNNNNNNUNNNUNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNN', 
                'iter_no=1000', 'seed=CUGGGGUGUCUAGCAGAUUUCAUAUGUCUGGUAUCAAUUCUGUGUUAGCUUUAAGCGCAAUCGCCCCAG']
        '''
        logging.error("Usage: python rna_designer.py <target_structure> <target_sequence> [energy=<target_energy>] "
                      "[mr=<target_mr>] [seed=<seed sequence>] [iter_no=<No iteretions>")
        sys.exit(-1)
        '''
    params = {'structure': argv[1], 'sequence': argv[2]}
    extended_args = argv[3:]
    for arg in extended_args:
        extended = arg.split('=')
        params[extended[0].strip()] = extended[1].strip()
    print(run_rnafbinv(**params))
