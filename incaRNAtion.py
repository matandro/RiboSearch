from tempfile import NamedTemporaryFile as NTF
from subprocess import Popen, PIPE
import os
import re
import logging
import sys

#INCARNATION_PATH = "/opt/algorithm/incaRNAtion2/IncaRNAtion"
INCARNATION_PATH = "/opt/algorithm/incaRNAtion/IncaRNAtion-master/src/IncaRNAtion.py"
PYTHON_PATH = "python"
MAX_ATTEMPT = 5


def run_incaRNAtion(structure, amount_to_generate, gc_content = 0.5, sequence_constraints = None):
    result = []
    attempt = 0
    if sequence_constraints is None:
        sequence_constraints = 'N' * len(structure)
    tmp_file = NTF(dir='.', delete=False)
    try:
        tmp_file.write('{}\n'.format(structure).encode())
        tmp_file.close()
        #param_list = [INCARNATION_PATH, '-a', '1', '-d', tmp_file.name,
        param_list = [PYTHON_PATH, INCARNATION_PATH, '-a', '1', '-d', tmp_file.name,
                      '-c', sequence_constraints, '-no_profile', '-s_gc', str(gc_content)]
        while len(result) < amount_to_generate and attempt < MAX_ATTEMPT:
            attempt += 1
            param_list.append(str(amount_to_generate - len(result)))
            result.extend(_single_run(param_list))
            result = list(set(result)) # remove duplicates
            param_list.pop()
    finally:
        os.remove(tmp_file.name)
    return result


#RES_REGEXP = re.compile(r'(?P<sequence>[UGCA]+) \((?P<energy>.)*\)')
RES_REGEXP = re.compile(r'(?P<sequence>[UGCA]+)')


def _single_run(param_list):
    def analyze_res(res: str):
        full_result = []
        for line in res.split('\n'):
            match = RES_REGEXP.match(line)
            if match is not None:
                full_result.append(match.group('sequence'))
        return full_result
    try:
        with Popen(param_list, stdout=PIPE, stdin=PIPE) as proc:
            logging.debug("Running incaRNAtion on params: {}".format(param_list))
            result = proc.communicate()[0].decode()
    except:
        logging.error("Failed to run: '{}'".format(param_list))
        sys.exc_info()
        sys.exit(1)
    return analyze_res(result)


def generate_seeds(structure, amount_to_generate):
    return run_incaRNAtion(structure, amount_to_generate)


if __name__ == '__main__':
    print(generate_seeds('(((...(((...((((((...)))......(((...)))...(((...(((...)))...)))...)))...)))...)))', 10))
