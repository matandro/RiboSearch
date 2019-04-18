#!/usr/bin/env python3
"""
Generates statistics and images for dive results.
F
"""
from typing import List, Tuple, Optional, Dict, TextIO
import os
import shutil
from subprocess import Popen, PIPE
import logging
import infernal
from rnafbinv import vienna, shapiro_tree_aligner
from tempfile import NamedTemporaryFile as NTF


def analyze_sto(sto_file: str, name: str) -> Tuple[str, str, Dict[str, str]]:
    structure = ''
    sequence = ''
    tseq_align = ''
    pp_align = ''
    struct_align = ''
    qseq_align = ''
    sel_name = name
    with open(sto_file, 'r') as stk_input:
        for line in stk_input:
            name_line = '#=GR {}'.format(name)
            if line.startswith(sel_name):
                # sel_name adds support for file with multiple alignment (take first which is best)
                sel_name = line.split(maxsplit=1)[0].strip()
                tseq_align += line.split(maxsplit=1)[1].strip()
            elif line.startswith(name_line):
                pp_align += line.rsplit('PP', 1)[1].strip()
            elif line.startswith('#=GC SS_cons'):
                struct_align += line[12:].strip()
            elif line.startswith('#=GC RF'):
                qseq_align += line[7:].strip()
        origin_map_for = {}
        origin_map_rev = {}
        origin_new_map = {}
        new_origin_map = {}
        bracket_stack = []
        if pp_align == "":
            pp_align = '.' * len(tseq_align)
        index = 0
        for tseq_char, pp_char, struct_char, qseq_char in zip(tseq_align, pp_align, struct_align, qseq_align):
            if struct_char in '[{(<':
                bracket_stack.append(index)
            elif struct_char in ']})>':
                close_index = bracket_stack.pop()
                origin_map_for[close_index] = index
                origin_map_rev[index] = close_index
            if tseq_char.upper() in 'AGCU':
                if struct_char in '[]{}()<>':
                    origin_new_map[len(structure)] = index
                    new_origin_map[index] = len(structure)
                structure += struct_char
                sequence += tseq_char.upper()
            index += 1
        fix_structure = ''
        for index, struct_char in enumerate(structure):
            if struct_char in '[{(<':
                open_index = origin_new_map.get(index)
                close_index = origin_map_for.get(open_index)
                if new_origin_map.get(close_index) is None:
                    fix_structure += '.'
                else:
                    fix_structure += '('
            elif struct_char in ']})>':
                close_index = origin_new_map.get(index)
                open_index = origin_map_rev.get(close_index)
                if new_origin_map.get(open_index) is None:
                    fix_structure += '.'
                else:
                    fix_structure += ')'
            else:
                fix_structure += '.'
        structure = fix_structure
    align = {'tseq': tseq_align, 'pp': pp_align, 'struct': struct_align, 'qseq': qseq_align}
    return structure, sequence, align


# CMSEARCH for aligned structure (assuming seqdb_path is a fasta file with a single sequence
def get_cm_struct(cm_file_path: str, seqdb_path: str, debug: bool=False) \
        -> Tuple[str, str]:
    structure = None
    sequence = None
    temp_out = None
    try:
        temp_out = NTF(dir='.', delete=False)
        temp_out.close()
        param_list = [os.path.join(infernal.INFENRAL_PATH, infernal.CMSEARCH_EXE), '--max', '--incE', '10', '-g',
                      '-A', temp_out.name, cm_file_path, seqdb_path]
        logging.info("Starting cm search {}".format(param_list))
        with Popen(param_list, stdout=PIPE, stdin=PIPE) as proc:
            output, err = proc.communicate()
            ret_code = proc.wait()
            if ret_code < 0:
                raise Exception(err)
            name = None
            all_file = ""
            with open(temp_out.name, 'r') as test_file:
                for line in test_file:
                    all_file += line
                    if len(line.strip()) > 1 and line[0] != '#' and name is None:
                        name = line.split()[0]
                        break
            if name is None:
                raise Exception("failed to figure new name\n{}".format(all_file))
            structure, sequence, align = analyze_sto(temp_out.name, name)
    except Exception as e:
        logging.error("Failed to search cm file {} on sequence db {}. ERROR: {}"
                      .format(cm_file_path, seqdb_path, e))
    finally:
        if temp_out is not None and os.path.exists(temp_out.name):
            if not debug:
                os.remove(temp_out.name)
            else:
                logging.info("Finished debug run on cm: {} output: {}".format(cm_file_path, temp_out.name))
    return structure, sequence


def generate_r2r(folder_path: str, sto_file_name: str, force: bool = False):
    pre_comp_sto = os.path.join(folder_path, "{}_cons.sto".format(sto_file_name[:-4]))
    meta_file_path = os.path.join(folder_path, "{}.meta".format(sto_file_name[:-4]))
    result_img_path = os.path.join(folder_path, "{}.svg".format(sto_file_name[:-4]))
    if os.path.exists(result_img_path) and not force:
        return result_img_path
    try:
        args = ["r2r", "--GSC-weighted-consensus", os.path.join(folder_path, sto_file_name), pre_comp_sto,
                "3", "0.97", "0.9", "0.75", "4", "0.97", "0.9", "0.75", "0.5", "0.1"]
        with Popen(args, stdout=PIPE, stdin=PIPE) as proc:
            proc.communicate()
        if not os.path.exists(pre_comp_sto) or os.stat(pre_comp_sto).st_size == 0:
            logging.error("Failed to compute consensus structure for {}".format(sto_file_name))
            return None
        data = None
        with open(pre_comp_sto, 'r') as sto_file:
            data = sto_file.readlines()
        with open(pre_comp_sto, 'w') as sto_file:
            for i in range(0, len(data))[::-1]:
                if data[i].startswith('//'):
                    data.insert(i, '#=GF R2R SetDrawingParam autoBreakPairs true\n')
            sto_file.writelines(data)
        with open(meta_file_path, 'w') as meta_file:
            meta_file.write(pre_comp_sto)
            meta_file.flush()
        args = ["r2r", "--disable-usage-warning", meta_file_path, result_img_path]
        with Popen(args) as proc:
            proc.communicate()
        if not os.path.exists(result_img_path) or os.stat(result_img_path).st_size == 0:
            if os.path.exists(result_img_path) and os.stat(result_img_path).st_size == 0:
                os.remove(result_img_path)
            logging.error("Failed to generate svg image for {}".format(result_img_path))
            return None
    finally:
        if os.path.exists(pre_comp_sto):
            os.remove(pre_comp_sto)
        if os.path.exists(meta_file_path):
            os.remove(meta_file_path)
    logging.info("Created r2r for {}: {}".format(sto_file_name, result_img_path))
    return result_img_path


def generate_cmv(folder_path: str, cm_file_name: str, sto_file_name: str,
                 model_type: str = 'detailed', force: bool = False):
    out_svg = os.path.join(folder_path, "{}_{}.svg".format(cm_file_name[:-3], model_type))
    if os.path.exists(out_svg) and not force:
        return out_svg
    args = ["conda", "run", "CMV", "-m", os.path.join(folder_path, cm_file_name),
            '-s', os.path.join(folder_path, sto_file_name), '-d', model_type, '-f', 'svg',
            '-p', folder_path]
    with Popen(args) as proc:
        proc.communicate()
    normal_out_svg = os.path.join(folder_path, "{}.svg".format(cm_file_name[6:-3], model_type))
    if os.path.exists(normal_out_svg):
        shutil.move(normal_out_svg, out_svg)
    else:
        return None
    return out_svg


def generate_visuals(folder_path: str):
    files = os.listdir(folder_path)
    cm_files = [file_name for file_name in files if file_name[-3:] == '.cm']
    for cm_name in cm_files:
        sto_name = "{}.sto".format(cm_name[:-3])
        generate_cmv(folder_path, cm_name, sto_name)
        generate_r2r(folder_path, sto_name)


def generate_filter(folder_path: str, cuttoff: int) -> List[str]:
    result = []
    with open(os.path.join(folder_path, "FINAL_summary"), 'r') as sum_file:
        sum_file.readline()
        for line in sum_file:
            values = line.strip().split('\t')
            if len(values) < 7:
                continue
            if int(values[2]) > cuttoff:
                result.append(values[0].strip())
    with open(os.path.join(folder_path, "Filter_list_{}.txt".format(cuttoff)), 'w') as out_file:
        out_file.write("\n".join(result))
        out_file.flush()
    return result


def generate_common(folder_path):
    def add_count(main_map: Dict[str, Dict[str, int]], key: str, value: str):
        count_map = main_map.get(key, {})
        count = count_map.get(value, 0)
        count_map[value] = count + 1
        main_map[key] = count_map

    def print_counts(main_map: Dict[str, Dict[str, int]], out_file: TextIO):
        for first_col, count_map in main_map.items():
            for second_col, count in count_map.items():
                out_file.write('{}\t{}\t{}\n'.format(first_col, second_col, count))

    sequence_map = {}
    location_map = {}
    with open(os.path.join(folder_path, 'FINAL_all'), 'r') as all_file, \
            open(os.path.join(folder_path, 'sequence_rep'), 'w') as sequence_rep_file, \
            open(os.path.join(folder_path, 'location_rep'), 'w') as location_rep_file:
        all_file.readline()
        sequence_rep_file.write('sequence\tdesign_code\tcount\n')
        location_rep_file.write('identifier\tdesign_code\tcount\n')
        for line in all_file:
            # 0:design_code, 1:identifier, 2:score, 3:E-value, 4:sequence
            items = line.strip().split('\t')
            design_code = items[0].strip()
            identifier = items[1].strip()
            sequence = items[4].strip()
            add_count(sequence_map, sequence, design_code)
            add_count(location_map, identifier, design_code)
        print_counts(sequence_map, sequence_rep_file)
        print_counts(location_map, location_rep_file)
        sequence_rep_file.flush()
        location_rep_file.flush()


def get_structures(folder_path: str, sequence: str, identifier: str, design_code: str) -> Tuple[str, str, str, str]:
    seq_db = None
    cm_struct = None
    min_struct = None
    cm_min_struct = None
    new_sequence = None
    cm_path = os.path.join(folder_path, "FINAL_{}.cm".format(design_code))
    try:
        seq_db = infernal.generate_fasta({identifier: sequence})
        cm_struct, new_sequence = get_cm_struct(cm_path, seq_db.name)
        if len(sequence) != len(new_sequence):
            logging.error("Sequence structure length differ\n{}\n{}\n{}".format(sequence, new_sequence, cm_struct))
        min_struct = vienna.fold(new_sequence)['MFE']
        cm_min_struct = vienna.fold(new_sequence, structure_constraints=cm_struct)['MFE']
    except Exception as e:
        print("Error on ID: {}\n{}".format(identifier, e))
    finally:
        if seq_db is not None and os.path.exists(seq_db.name):
            os.remove(seq_db.name)
    return new_sequence, cm_struct, min_struct, cm_min_struct


def extend_final_all(folder_path: str):
    def search_cm(target_seq: str, target_ident: str) -> Tuple[bool, Optional[str]]:
        seq_db = None
        try:
            seq_db = infernal.generate_fasta({target_ident: target_seq})
            cm_res = infernal.search_cm(cm_purine_path, seq_db.name, inc_e=10)
        finally:
            if seq_db is not None and os.path.exists(seq_db.name):
                os.remove(seq_db.name)
        if cm_res is not None and cm_res:
            return True, cm_res[0]['E-value']
        else:
            return False, None

    def calc_score(source_sequence: str, source_structure: str) -> float:
        source_tree = shapiro_tree_aligner.get_tree(source_structure, source_sequence)
        _, score = shapiro_tree_aligner.align_trees(source_tree, target_tree)
        return score

    cm_purine_path = os.path.join(folder_path, 'RF00167.cm')
    target_tree = shapiro_tree_aligner.get_tree("((((((((...(.(((((.......))))).)........((((((.......))))))..))))))))",
                                                "NNNNNNNNUNNNNNNNNNNNNNNNNNNNNNNNNUNNNUNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNN")
    with open(os.path.join(folder_path, 'FINAL_all'), 'r') as all_in_file, \
            open(os.path.join(folder_path, 'FINAL_all_ext'), 'w') as all_out_file:
        header = all_in_file.readline().strip()
        header += "\tis_rfam_matched\trfam_evalue\talignment_score_min\talignment_score_cm\talignment_score_cm_min\n"
        all_out_file.write(header)
        for line in all_in_file:
            # 0:design_code, 1:identifier, 2:score, 3:E-value, 4:sequence, 5:round
            items = line.strip().split('\t')
            design_code = items[0].strip()
            identifier = items[2].strip()
            sequence = items[4].strip()
            # check rfam purine matching
            is_found, evalue = search_cm(sequence, identifier)
            items.append('1' if is_found else '0')
            items.append(str(evalue))
            # get scores
            new_sequence, cm_struct, min_struct, cm_min_struct = get_structures(folder_path, sequence,
                                                                                identifier, design_code)
            items.append(str(calc_score(new_sequence, min_struct)))
            items.append(str(calc_score(new_sequence, cm_struct)))
            items.append(str(calc_score(new_sequence, cm_min_struct)))
            all_out_file.write('{}\n'.format('\t'.join(items)))


def generate_statistics(folder_path: str):
    generate_common(folder_path)
    extend_final_all(folder_path)


if __name__ == '__main__':
    threshold = [5, 10, 25, 50]
    base_folder = "/opt/home/matan/Dropbox/PHd/RiboSearch/IncManual"
    for cut in threshold:
        left_list = generate_filter(base_folder, cut)
    generate_visuals(base_folder)
    generate_statistics(base_folder)
