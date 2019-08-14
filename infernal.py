from subprocess import Popen, PIPE
import logging
import os
from rnafbinv import vienna
from tempfile import NamedTemporaryFile as NTF
from typing import List, Dict
from Bio import SeqIO
from enum import Enum


INFENRAL_PATH = "/opt/algorithm/infernal/bin/"
CMBUILD_EXE = "cmbuild"
CMSEARCH_EXE = "cmsearch"
CMCALIBRATE_EXE = "cmcalibrate"
CMALIGN_EXE = "cmalign"
SHELL_SEQ_SCRIPT = "/opt/algorithm/RiboSearch/infernal_pull_seq.sh"
ESL_FETCH = "/opt/algorithm/infernal/bin/esl-sfetch"
#SHELL_SEQ_SCRIPT = "infernal_pull_seq.sh"
STOCKHOLM_FORMAT = "# STOCKHOLM 1.0"
FASTA_LINE_LENGTH = 80


def generate_fasta(sequences: Dict[str, str]) -> NTF:
    tmp_file = NTF(mode='w+', dir='.', delete=False, encoding="utf-8")
    for topic, sequence in sequences.items():
        tmp_file.write('> {}\n'.format(topic))
        for fasta_line in [sequence[i:i+FASTA_LINE_LENGTH] for i in range(0, len(sequence), FASTA_LINE_LENGTH)]:
            tmp_file.write('{}\n'.format(fasta_line))
    tmp_file.close()
    return tmp_file


def align_sequences(sequences: Dict[str, str], cm_path: str, out_align_path: str, in_align_path: str=None,
                    timeout: int=None) -> bool:
    result = False
    tmp_sequences_fasta = generate_fasta(sequences)
    try:
        param_list = [os.path.join(INFENRAL_PATH, CMALIGN_EXE), '-g', '-o', out_align_path]
        if in_align_path is not None:
            param_list += ['--mapali', in_align_path]
        param_list += [cm_path, tmp_sequences_fasta.name]
        logging.info("Aligning sequences to CM file {}".format(cm_path))
        with Popen(param_list, stdout=PIPE, stdin=PIPE) as proc:
            _, stderr_out = proc.communicate(timeout=timeout)
            if stderr_out is not None:
                logging.error(stderr_out)
            exit_code = proc.wait()
        if os.path.exists(out_align_path) and exit_code == 0:
            result = True
    except Exception as e:
        if os.path.exists(out_align_path):
            os.remove(out_align_path)
        logging.error("Failed to sequence to cm file {}. ERROR: {}".format(cm_path, e))
    finally:
        if tmp_sequences_fasta is not None and os.path.exists(tmp_sequences_fasta.name):
            os.remove(tmp_sequences_fasta.name)
    return result


def generate_stockholm(sequence: str, structure: str=None) -> NTF:
    if structure is None:
        structure = vienna.fold(sequence)['MFE']
    tmp_file = NTF(mode='w+', dir='.', delete=False, encoding="utf-8")
    tmp_file.write('{}\n'.format(STOCKHOLM_FORMAT))
    tmp_file.write('seq1\t{}\n'.format(sequence))
    tmp_file.write('#=GC SS_cons\t{}\n//'.format(structure))
    tmp_file.close()
    return tmp_file


def generate_cm(stockholm_path: str, outcm_path: str, cpus: int=None, force: bool=False) -> bool:
    result = False
    try:
        if not os.path.exists(outcm_path) or force:
            param_list = [os.path.join(INFENRAL_PATH, CMBUILD_EXE), '-F']
            if cpus is not None:
                param_list += ['--cpu', str(cpus)]
            param_list += [outcm_path, stockholm_path]
            logging.info("Generating CM file {} for file: {}".format(outcm_path, stockholm_path))
            with Popen(param_list, stdout=PIPE, stdin=PIPE) as proc:
                ret_code = proc.wait()
                if ret_code < 0:
                    raise Exception("cmbuild ended with error code {}".format(ret_code))
        if not is_calibrated(outcm_path):
            param_list = [os.path.join(INFENRAL_PATH, CMCALIBRATE_EXE)]
            if cpus is not None:
                param_list += ['--cpu', str(cpus)]
            param_list.append(outcm_path)
            logging.info("Calibrating CM file {}".format(outcm_path))
            with Popen(param_list, stdout=PIPE, stdin=PIPE) as proc:
                proc.communicate()
                ret_code = proc.wait()
                if ret_code < 0:
                    raise Exception("cmcalibrate ended with error code {}".format(ret_code))
        result = os.path.exists(outcm_path)
        logging.info("Finished CM file")
    except Exception as e:
        if os.path.exists(outcm_path):
            os.remove(outcm_path)
        logging.error("Failed to generate and calibrate cm file {}. ERROR: ".format(outcm_path, e))
        raise e
    return result


def generate_single_seq_cm(sequence: str, outcm_path: str, structure: str=None, cpus: int=None) -> bool:
    result = False
    tmp_stockholm = generate_stockholm(sequence, structure)
    try:
        result = generate_cm(tmp_stockholm.name, outcm_path, cpus=cpus)
    finally:
        if tmp_stockholm is not None and os.path.exists(tmp_stockholm.name):
            os.remove(tmp_stockholm.name)
    return result


def fetch_seq_tlbout(tlbout_path: str, fasta_orig: str) -> List[Dict[str, str]]:
    results = []
    ssi_file = None
    try:
        ssi_file = "{}.ssi".format(fasta_orig)
        if os.path.exists(ssi_file):
            # exists before then don't create and remove
            ssi_file = None
        else:
            param_list = [ESL_FETCH, '--index', fasta_orig]
            logging.info("generating index {}".format(param_list))
            with Popen(param_list, stdout=PIPE, stdin=PIPE) as proc:
                proc.communicate()
        temp_fasta_out = NTF(dir='.', delete=False)
        temp_fasta_out.close()
        param_list = ['sh', SHELL_SEQ_SCRIPT, tlbout_path, fasta_orig, temp_fasta_out.name]
        logging.info("Retrieving sequence {}".format(param_list))
        with Popen(param_list, stdout=PIPE, stdin=PIPE) as proc:
            proc.communicate()
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
        if ssi_file is not None and os.path.exists(ssi_file):
            os.remove(ssi_file)
    return results


def output_search_analyze_tlbout(output: str) -> List[Dict[str, str]]:
    search_list = []
    for line in output.split('\n'):
        line = line.strip()
        if line[:1] != "#" and len(line) > 0:
            words = line.split()
            single_res = {'target name': words[0], 'model from': words[5], 'model to': words[6],
                          'seq from': words[7], 'seq to': words[8], 'strand': words[9],
                          'score': words[14], 'E-value': words[15], 'description': words[17],
                          'identifier': '{}/{}-{}'.format(words[0], words[7], words[8])}
            search_list.append(single_res)
    return search_list


def output_search_analyze(output: str) -> List[Dict[str, str]]:
    return output_search_analyze_tlbout(output)
    # return output_search_analyze_normal(output)


# NOT USED
def output_search_analyze_normal(output: str) -> List[Dict[str, str]]:
    def split_results(hit_section_text: str) -> List[str]:
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
            res_seq = res_seq.replace('-', '').upper()
            res = {'target name': res_name, 'score': res_info[3], 'E-value': res_info[2], 
                   'seq from': res_info[9], 'seq to': res_info[10], 'strand': res_info[11],
                   'model from': res_info[6], 'model to': res_info[7], 'sequence': res_seq,
                   'identifier': '{}/{}-{}'.format(res_name, res_info[9], res_info[10])}
            search_list.append(res)
    return search_list


def is_calibrated(cm_file_path: str) -> bool:
    calibrated = False
    with open(cm_file_path, 'r') as cm_file:
        for line in cm_file:
            if line[0:3] == 'COM' and 'cmcalibrate' in line:
                calibrated = True
                break
    return calibrated


class ResType(Enum):
    ERIC = 1
    TBLOUT = 2
    MANUAL = 3


def search_cm(cm_file_path: str, seqdb_path: str, debug: bool=False,
              res_type: ResType=ResType.ERIC, inc_e: float=None, cpus: int=None) -> List[Dict[str, str]]:
    def merge_eric(table_results: List[Dict[str, str]], eric_results: List[Dict[str, str]]):
        for eric_res, table_res in zip(eric_results, table_results):
            eric_target = eric_res.get("target name")
            target, loc_str = eric_target.rsplit('/', 1)
            seq_from, seq_to = loc_str.split('-', 1)
            if table_res.get('target name') != target or table_res.get('seq from') != seq_from \
                    or table_res.get('seq to') != seq_to:
                logging.error("Rows do not match: table - {} eric - {}".format(table_res, eric_res))
            if table_res.get('target name') == target:
                table_res['sequence'] = eric_res.get('sequence')

    results = None
    temp_out = None
    try:
        temp_out = NTF(dir='.', delete=False)
        temp_out.close()
        param_list = [os.path.join(INFENRAL_PATH, CMSEARCH_EXE), '--tblout', temp_out.name]
        if inc_e is not None:
            param_list += ['--incE', str(inc_e)]
        if cpus is not None:
            param_list += ['--cpu', str(cpus)]
        param_list += [cm_file_path, seqdb_path]
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
                    logging.warning("Something strange in infernal result analysis. tlbout lines = {},"
                                    " esl_sfetch lines = {}\n{}\n{}".format(len(tlbout_results), len(results),
                                                                            tlbout_results, results))
                else:
                    merge_eric(tlbout_results, results)
                    results = tlbout_results
            elif res_type == ResType.TBLOUT:
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
    test_sequence = "AGGUCCUAGUGCAGCGGGACUUUUUUUCUAAAGUCGUUGAGAGGAGGAGUCGUCAGACCAGAUAGCUUUGAUGUCCUGAU" \
                    "CGGAAGGAUCGUUGGCCCCC"
    test_sequence = "GGAGGCCGCUUGCCCUCC"
    is_build_cm = generate_single_seq_cm(test_sequence, "test.cm")
    print("is build CM {}".format(is_build_cm))
    if is_build_cm:
        print(search_cm("test.cm", "test_seq_db.fasta", debug=True))
        # print(search_cm("test.cm", "no_find_seq_db.fasta", debug=True))
        # print(search_cm("1_3.cm", "/DB/fasta_db/Tbrucei/Tbrucei", debug=True))
