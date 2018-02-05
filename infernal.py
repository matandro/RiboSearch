from subprocess import Popen, PIPE
import logging
import os
import vienna
from tempfile import NamedTemporaryFile as NTF
import re
from Bio import SeqIO
from enum import Enum


INFENRAL_PATH = "/opt/algorithm/infernal/bin/"
CMBUILD_EXE = "cmbuild"
CMSEARCH_EXE = "cmsearch"
CMCALIBRATE_EXE = "cmcalibrate"
SHELL_SEQ_SCRIPT = "/opt/algorithm/RiboSearch/infernal_pull_seq.sh"
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


def fetch_seq_tlbout(tlbout_path, fasta_orig):
    results = []
    try:
        temp_fasta_out = NTF(dir='.', delete=False)
        temp_fasta_out.close()
        param_list = ['sh', SHELL_SEQ_SCRIPT, tlbout_path, fasta_orig, temp_fasta_out.name]
        logging.info("Retrieving sequence {}".format(param_list))
        with Popen(param_list, stdout=PIPE, stdin=PIPE) as proc:
            proc.wait()
        if os.path.exists(temp_fasta_out.name):
            for seq_record in SeqIO.parse(temp_fasta_out.name, "fasta"):
                results.append({"target name": seq_record.id, "sequence": "{}".format(seq_record.seq)})
        else:
            raise Exception("fasta file not created")
    except Exception as exc:
        logging.error("Failed to retrieve sequences {}".format(exc))
    finally:
        if temp_fasta_out is not None and os.path.exists(temp_fasta_out.name):
            os.remove(temp_fasta_out.name)
    return results


def output_search_analyze_tlbout(output):
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


def output_search_analyze(output):
    return output_search_analyze_tlbout(output)
    #return output_search_analyze_normal(output)

# NOT USED
def output_search_analyze_normal(output):
    def split_results(hit_section_text):
        res = []
        collect = None
        for line in hit_section_text.split('\n'):
            if len(line) > 0 and line[0] == '>':
                if collect is not None:
                    res.append(collect)
                collect = ''
            if collect is not None:
                if not collect == '':
                    collect += '\n'
                collect += line
        if collect is not None:
            res.append(collect)
        return res
    search_list = []
    start_index = output.find('Hit alignments:')
    end_index = output.find('Internal CM pipeline statistics summary:')
    output = output[start_index:end_index].split('\n', 1)[1]
    if 'No hits detected that satisfy reporting thresholds' not in output:
        result_list = split_results(output)
        for single_hit in result_list:
            res_name = single_hit.split('\n',1)[0].strip('>')
            res_section= single_hit.split('\n\n')
            res_info = [ainfo.strip() for ainfo in res_section[0].split('\n')[3].split()]
            res_section = res_section[1:]
            res_seq = ''
            for partial_info in res_section:
                lines = partial_info.split('\n')
                index_nc = -1
                for i in range(0,len(lines)):
                    if lines[i].strip() == 'NC':
                        index_nc = i
                        break
                if index_nc >= 0:
                    temp_split = partial_info.split('\n')[index_nc + 4].rsplit(' ',1)[0].strip()
                    res_seq += temp_split.split(' ',3)[2]
            res_seq = res_seq.replace('-','').upper()
            res = {'target name': res_name, 'score': res_info[3], 'E-value': res_info[2], 
                   'seq from': res_info[9], 'seq to': res_info[10], 'strand': res_info[11],
                   'model from': res_info[6], 'model to': res_info[7], 'sequence': res_seq}
            search_list.append(res)
    return search_list


def is_calibrated(cm_file_path):
    calibrated = False
    with open(cm_file_path, 'r') as cm_file:
        for line in cm_file:
            if line[0:3]  == 'COM' and 'cmcalibrate' in line:
                calibrated = True
                break
    return calibrated


class ResType(Enum):
    ERIC = 1
    TLBOUT = 2
    MANUAL = 3


def search_cm(cm_file_path, seqdb_path, debug=False, res_type=ResType.ERIC):
    results = None
    temp_out = None
    try:
        temp_out = NTF(dir='.', delete=False)
        temp_out.close()
        param_list = [os.path.join(INFENRAL_PATH, CMSEARCH_EXE), '--tblout', #'-A', 'something.stk', '-o',
                      temp_out.name, cm_file_path, seqdb_path]
        logging.info("Starting cm search {}".format(param_list))
        with Popen(param_list, stdout=PIPE, stdin=PIPE) as proc:
            output, err = proc.communicate()
            ret_code = proc.wait()
            if ret_code < 0:
                raise Exception(err)
            # keeping this to compare for errors!
            with open(temp_out.name, 'r') as output_file:
                results = ''
                for line in output_file:
                    results += line
            tlbout_results = output_search_analyze(results)
            if res_type == ResType.ERIC:
                # getting actual results (sequnces)
                results = fetch_seq_tlbout(temp_out.name, seqdb_path)
                if len(tlbout_results) != len(results):
                    logging.warning("Something strange in infernal result analysis. tlbout lines = {}, esl_sfetch lines = {}".format(len(tlbout_results), len(results)))
            elif res_type == ResType.TLBOUT:
                results = tlbout_results
            logging.info("Finisied cm search {} results".format(len(results)))
    except Exception as e:
        logging.error("Failed to search cm file {} on sequence db {}. ERROR: {}"
                      .format(cm_file_path, seqdb_path, e))
    finally:
        if temp_out is not None and os.path.exists(temp_out.name):
            if not debug:
                os.remove(temp_out.name)
            else:
                logging.info("Finished debug run on cm: {} output: {}".format(cm_file_path, temp_out.name))
    return results


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_sequence = "AGGUCCUAGUGCAGCGGGACUUUUUUUCUAAAGUCGUUGAGAGGAGGAGUCGUCAGACCAGAUAGCUUUGAUGUCCUGAUCGGAAGGAUCGUUGGCCCCC" #"GGAGGCCGCUUGCCCUCC"
    is_build_cm = generate_cm(test_sequence,"test.cm")
    print("is build CM {}".format(is_build_cm))
    if is_build_cm:
#        print(search_cm("test.cm", "test_seq_db.fasta", debug=True))
#        print(search_cm("test.cm", "no_find_seq_db.fasta", debug=True))
        print(search_cm("1_3.cm", "/DB/fasta_db/Tbrucei/Tbrucei", debug=True))

