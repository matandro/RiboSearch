from subprocess import Popen, PIPE
import logging
import os
import vienna
from tempfile import NamedTemporaryFile as NTF


INFENRAL_PATH = "/opt/algorithm/infernal/bin/"
CMBUILD_EXE = "cmbuild"
CMSEARCH_EXE = "cmsearch"
CMCALIBRATE_EXE = "cmcalibrate"
STOCKHOLM_FORMAT = "# STOCKHOLM 1.0"
MAX_THREADS = 4


def generate_stockholm(sequence):
    structure = vienna.fold(sequence)['MFE']
    tmp_file = NTF(dir='.', delete=False)
    tmp_file.write('{}\n'.format(STOCKHOLM_FORMAT).encode())
    tmp_file.write('seq1\t{}\n'.format(sequence).encode())
    tmp_file.write('#=GC SS_cons\t{}\n//'.format(structure).encode())
    tmp_file.close()
    return tmp_file


def generate_cm(sequence, outcm_path):
    result = False
    tmp_stockholm = generate_stockholm(sequence)
    try:
        if not os.path.exists(outcm_path):
            param_list = [os.path.join(INFENRAL_PATH, CMBUILD_EXE), '-F', outcm_path,
                          tmp_stockholm.name]
            logging.info("Generating CM file {} for sequence: {}".format(outcm_path, sequence))
            with Popen(param_list, stdout=PIPE, stdin=PIPE) as proc:
                ret_code = proc.wait()
                if ret_code < 0:
                    raise Exception()
            param_list = [os.path.join(INFENRAL_PATH, CMCALIBRATE_EXE), '--cpu', str(MAX_THREADS), outcm_path]
        if not is_calibrated(outcm_path):
            logging.info("Calibrating CM file")
            with Popen(param_list, stdout=PIPE, stdin=PIPE) as proc:
                ret_code = proc.wait()
                if ret_code < 0:
                    raise Exception()
            result = os.path.exists(outcm_path)
        logging.info("Finished CM file")
    except Exception as e:
        if os.path.exists(outcm_path):
            os.remove(outcm_path)
        logging.error("Failed to generate and calibrate cm file {}. ERROR: ".format(outcm_path, e))
    finally:
        if tmp_stockholm is not None and os.path.exists(tmp_stockholm.name):
            os.remove(tmp_stockholm.name)
    return result


def output_search_analyze(output):
    search_list = []
    for line in output.split('\n'):
        line = line.strip()
        if line[:1] != "#" and len(line) > 0:
            words = line.split()
            single_res = {'target name': words[0], 'model from': words[5], 'model to': words[6],
                          'seq from': words[7], 'seq to': words[8], 'strand': words[9],
                          'score': words[14], 'E-value': words[15], 'description': words[17]}
            search_list.append(single_res)
    return search_list


def is_calibrated(cm_file_path):
    calibrated = False
    with open(cm_file_path, 'r') as cm_file:
        for line in cm_file:
            if line[0:3]  == 'COM' and 'cmcalibrate' in line:
                calibrated = True
                break
    return calibrated



def search_cm(cm_file_path, seqdb_path):
    results = None
    temp_out = None
    try:
        temp_out = NTF(dir='.', delete=False)
        temp_out.close()
        param_list = [os.path.join(INFENRAL_PATH, CMSEARCH_EXE), '--tblout', temp_out.name,
                      cm_file_path, seqdb_path]
        logging.info("Starting cm search {}".format(param_list))
        with Popen(param_list, stdout=PIPE, stdin=PIPE) as proc:
            output, err = proc.communicate()
            ret_code = proc.wait()
            if ret_code < 0:
                raise Exception(err)
            with open(temp_out.name, 'r') as output_file:
                results = ''
                for line in output_file:
                    results += line
            results = output_search_analyze(results)
            logging.info("Finisied cm search {} results".format(len(results)))
    except Exception as e:
        logging.error("Failed to search cm file {} on sequence db {}. ERROR: {}"
                      .format(cm_file_path, seqdb_path, e))
    finally:
        if temp_out is not None and os.path.exists(temp_out.name):
            os.remove(temp_out.name)
    return results


if __name__ == "__main__":
    test_sequence = "GGAGGCCGCUUGCCCUCC"
    rebuild = False
    if rebuild:
        is_build_cm = generate_cm(test_sequence,"test.cm")
        print("is build CM {}".format(is_build_cm))
    if not rebuild or is_build_cm:
        print(search_cm("test.cm", "test_seq_db.fasta"))
