#!/usr/bin/env python3
"""
Generates statistics and images for dive results.
F
"""
from collections import OrderedDict
from typing import List, Tuple, Dict, TextIO, Optional, Set
import compare_rfam
import os
import shutil
from subprocess import Popen, PIPE
import logging
import infernal
from rnafbinv import vienna, shapiro_tree_aligner
from tempfile import NamedTemporaryFile as NTF
from subprocess import Popen
import time
import result_dive
import gffutils
import filecmp


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


def generate_r2r(folder_path: str, sto_file_name: str, force: bool = False, output_name: Optional[str]= None):
    pre_comp_sto = os.path.join(folder_path, "{}_cons.sto".format(sto_file_name[:-4]))
    meta_file_path = os.path.join(folder_path, "{}.meta".format(sto_file_name[:-4]))
    result_img_path = os.path.join(folder_path, "{}.r2r.svg".format(sto_file_name[:-4]))
    if output_name is not None:
        result_img_path = output_name
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
    args = ["/opt/home/matan/anaconda2/bin/conda", "run", "CMV", "-m", os.path.join(folder_path, cm_file_name),
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


def get_sto_targets(sto_path: str) -> Set[str]:
    results = set()
    with open(sto_path, 'r') as input_sto:
        for line in input_sto:
            stripped_line = line.strip()
            if stripped_line != "" and stripped_line[0] != '#':
                results.add(stripped_line.split(maxsplit=1)[0])
    return results


def recreate_cm(folder_path: str):
    # read inputs
    gather_designs = {}
    with open(os.path.join(folder_path, "FINAL_summary"), "r") as input_sum:
        input_sum.readline()
        for line in input_sum:
            if line.strip() == '':
                continue
            parts = line.strip().split('\t')
            gather_designs[parts[0]] = parts[5]
    gather_results = {}
    with open(os.path.join(folder_path, "FINAL_all"), 'r') as input_all:
        input_all.readline()
        for line in input_all:
            if line.strip() == '':
                continue
            parts = line.strip().split('\t')
            res_map = gather_results.get(parts[0], {})
            res_map[parts[1]] = parts[4]
            gather_results[parts[0]] = res_map
    # start calculations
    folder = vienna.LiveRNAfold()
    folder.start()
    for design_code, sequence in gather_designs.items():
        structure = folder.fold(sequence)['MFE']
        cm_path = os.path.join(folder_path, "{}.cm".format(design_code))
        sto_path = os.path.join(folder_path, "{}.sto".format(design_code))
        if os.path.exists(cm_path):
            continue
        temp_cm_path = "{}_tmp".format(cm_path)
        temp_sto_path = "{}_tmp".format(sto_path)
        if not infernal.generate_single_seq_cm(sequence, cm_path, structure):
            print("Could not generate single cm for {}".format(design_code))
            exit(-1)
        if not infernal.align_sequences({'{}'.format(design_code): sequence}, cm_path, sto_path):
            print("Could not generate single sto for {}".format(design_code))
            exit(-1)
        design_results = gather_results.get(design_code)
        no_found = 0
        temp_fasta = infernal.generate_fasta(design_results)
        while no_found < len(design_results):
            results = infernal.search_cm(cm_path, temp_fasta.name, inc_e= 10.0)
            sto_parts = {}
            sto_target = get_sto_targets(sto_path)
            for item in results:
                if item['target name'] not in sto_target:
                    sto_parts[item['target name']] = item['sequence']
            if len(sto_parts) == 0:
                print("ERROR: no new sequences found for {} maxed at {} sequences out of {} original\nListing: {}"
                      .format(design_code, len(sto_target), len(design_results), [res for res in design_results.keys()
                                                                                  if res not in
                                                                                  get_sto_targets(sto_path)]))
                break
            if not infernal.align_sequences(sto_parts, cm_path, temp_sto_path, in_align_path=sto_path):
                print("Could not generate sto for {}".format(design_code))
                exit(-1)
            if filecmp.cmp(sto_path, temp_sto_path, shallow=False):
                print("ERROR: {} missing codes: {}".format(design_code, [res for res in design_results.keys() if
                                                                         res not in get_sto_targets(sto_path)]))
                shutil.move(temp_sto_path, sto_path)
                break
            shutil.move(temp_sto_path, sto_path)
            if not infernal.generate_cm(sto_path, temp_cm_path):
                print("Could not generate cm for {}".format(design_code))
                exit(-1)
            shutil.move(temp_cm_path, cm_path)
            no_found = len(results)
        os.remove(temp_fasta.name)


def generate_visuals(folder_path: str):
    #recreate_cm(folder_path)
    files = os.listdir(folder_path)
    cm_files = [file_name for file_name in files if file_name[-3:] == '.cm']
    for cm_name in cm_files:
        sto_name = "{}.sto".format(cm_name[:-3])
        #generate_cmv(folder_path, cm_name, sto_name)
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
    cm_path = os.path.join(folder_path, "{}.cm".format(design_code))
    #cm_path = os.path.join(folder_path, "FINAL_{}.cm".format(design_code))
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


def get_gff(organism_id: str, BASE_DIR: str= '/opt/home/matan/Dropbox/PHd/RiboSearch/gff/', max_repeats: int=5,
            sleep_interval: int=10) -> Optional[str]:
    gff_path = os.path.join(BASE_DIR, "{}.gff".format(organism_id))
    DB = "nuccore"
    worked = False
    if not os.path.exists(gff_path):
        repeat = 0
        while repeat < max_repeats:
            with Popen(['wget', '-O', gff_path, 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?'
                                               'db={}&report=gff3&id={}'.format(DB, organism_id)]) as proc:
                proc.communicate()
                exit_code = proc.wait()
                if exit_code == 8:
                    if DB == 'nuccore':
                        DB = 'assembly'
                    else:
                        return None
                elif exit_code != 0:
                    repeat += 1
                else:
                    worked = True
            if not worked:
                time.sleep(sleep_interval)
            else:
                break
    else:
        worked = True
    if not worked:
        try:
            os.remove(gff_path)
        except:
            pass
        gff_path = None
    return gff_path


def tax_from_gff(gff_path: str) -> Optional[int]:
    result = None
    with open(gff_path, 'r') as input_gff:
        for line in input_gff:
            if line[0] != '#':
                break
            elif line.startswith('##species'):
                try:
                    result = int(line.split('id=')[1])
                except:
                    print("ERROR: failed to read taxonomy from species gff path {}".format(line.strip()))
    return result


def get_distance(start: int, end: int, target: gffutils.Feature) -> Tuple[int, str]:
    if start > end and target.strand == '+':
        raise ValueError("Mismatch of strands, comparing [{}-{}] to [{}-{}]({})".format(start, end, target.start,
                                                                                        target.end, target.strand))
    put_start = min(target.end, target.start)
    put_end = max(target.end, target.start)
    take_start = min(start, end)
    take_end = max(start, end)
    if take_end < put_start:
        distance = put_start - take_end
        relation = "near"
    elif take_start > put_end:
        distance = - (take_start - put_end)
        relation = "near"
    else:
        if take_start < put_start and take_end > put_end:
            #unlikely and sizes don't matter
            relation = "outside"
        elif take_start > put_start and take_end < put_end:
            from_end = min(take_start - put_start, put_end - take_end)
            relation = "inside {}{}".format('s' if from_end == (take_start - put_start) else 'e', from_end)
        elif take_end > put_start:
            relation = "partial e{}".format(take_end - put_start)
        elif take_start < put_start:
            relation = "partial s{}".format(put_start - take_start)
        distance = 0
    if target.strand == '-':
        distance = -distance
    return distance, relation


def change_dups(gff_path: str):
    line_map = OrderedDict()
    gff_out = '{}_{}'.format(gff_path, 'tmp')
    with open(gff_path, 'r') as gff_reader, \
        open(gff_out, 'w') as gff_write:
        for line in gff_reader:
            if line[0] == '#' or line.strip() == '':
                gff_write.write("{}\n".format(line.strip()))
                continue
            line_id = line.split('ID=', 1)[1].split(';', 1)[0]
            if line_map.get(line_id) is not None:
                num = 2
                while line_map.get("{}-{}".format(line_id, num)) is not None:
                    num += 1
                new_line_id = "{}-{}".format(line_id, num)
                new_line = line.replace("ID={}".format(line_id), "ID={}".format(new_line_id))
                line_map[new_line_id] = new_line
            else:
                line_map[line_id] = line
        for _, value in line_map.items():
            gff_write.write('{}\n'.format(value.strip()))
    shutil.move(gff_out, gff_path)


OUT_RANGE = 5000


def nearby_gene(gff_path: str, start: int, end: int, strand: str) -> Tuple[str, str, str]:
    gff_name = os.path.splitext(gff_path)[0]
    gff_db = "{}.db".format(gff_name)
    keep_building = True
    while keep_building:
        keep_building = False
        if not os.path.exists(gff_db):
            change_dups(gff_path)
            try:
                db = gffutils.create_db(gff_path, dbfn=gff_db, from_string=False)
            except ValueError as ve:
                if 'empty file provided' in str(ve):
                    try:
                        os.remove(gff_db)
                    except:
                        pass
                    return '', '', ''
                else:
                    print("failed to create {}".format(gff_db))
                    raise
            except:
                print("failed to create {}".format(gff_db))
                raise
        else:
            try:
                db = gffutils.FeatureDB(gff_db)
            except:
                os.remove(gff_db)
                keep_building = True
    gene_feature = None
    subname = ''
    product = ''
    for feature in db.region(start=max(0, min(start, end) - OUT_RANGE), end=max(start, end) + OUT_RANGE,
                             strand=strand):
        new_distance, new_relation = get_distance(start, end, feature)
        if feature.featuretype == 'region' and gene_feature is None:
            gene_feature = feature
            subname = '{} inside'.format(feature.featuretype)
            product = feature.attributes.get('mol_type')
        elif feature.featuretype == 'gene':
            if gene_feature is None or (gene_feature.featuretype == 'gene' and
                                        new_distance < get_distance(start, end, gene_feature)[0]):
                gene_feature = feature
                subname = '{} {} {}'.format(feature.featuretype, new_distance, new_relation)
                product = feature.attributes.get('product')
                if product is None:
                    product = feature.attributes.get('gene')
        elif feature.featuretype == 'CDS' and (gene_feature is None or gene_feature.featuretype == 'region' or
            feature.attributes.get('Parent') == gene_feature.attributes.get('ID')):
            f_product = feature.attributes.get('product')
            if f_product is not None:
                product = f_product
            # We have regions with CDS and no genes, otherwise gene has priority
            if gene_feature is not None and gene_feature.featuretype != 'gene':
                if gene_feature.featuretype == 'region' or \
                        (gene_feature.featuretype == 'CDS' and
                         new_distance < get_distance(start, end, gene_feature)[0]):
                    subname = '{} {} {}'.format(feature.featuretype, new_distance, new_relation)
                    gene_feature = feature
        elif feature.featuretype == 'five_prime_UTR' and new_distance == 0:
            subname = '{} {}'.format(feature.featuretype, new_relation)
            f_product = feature.attributes.get('product')
            if f_product is not None:
                product = f_product
        elif feature.featuretype == 'three_prime_UTR' and new_distance == 0:
            subname = '{} {}'.format(feature.featuretype, new_relation)
            f_product = feature.attributes.get('product')
            if f_product is not None:
                product = f_product
        elif feature.featuretype == 'exon' and new_distance == 0 and 'prime_UTR' not in subname:
            exon_num = ''
            fid_parts = feature.attributes.get('ID', default='')[0].rsplit('-', 1)
            if len(fid_parts) > 1:
                exon_num = fid_parts[1]
            subname = 'exon{} {}'.format(exon_num, new_relation)
            f_product = feature.attributes.get('product')
            if f_product is not None:
                product = f_product
    # get region tag if nothing was matched in in reverse
    if gene_feature is None:
        for feature in db.region(start=min(start, end), end=max(start, end)):
            if feature.featuretype == 'region':
                gene_feature = feature
                subname = '{} inside'.format(feature.featuretype)
                product = feature.attributes.get('mol_type')
                break
    gene_id = gene_feature.attributes.get('ID')
    return str(gene_id), str(product), subname


def extend_final_all(folder_path: str, continue_old: bool=True):
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
    skip_lines = 0
    old_exists = os.path.exists(os.path.join(folder_path, 'FINAL_all_ext'))
    if continue_old and old_exists:
        with open(os.path.join(folder_path, 'FINAL_all_ext'), 'r') as count_lines:
            count_lines.readline()
            for line in count_lines:
                skip_lines += 1
    elif old_exists:
        os.remove(os.path.join(folder_path, 'FINAL_all_ext'))
    with open(os.path.join(folder_path, 'FINAL_all'), 'r') as all_in_file, \
            open(os.path.join(folder_path, 'FINAL_all_ext'), 'a+') as all_out_file:
        header = all_in_file.readline().strip()
        header += "\tis_rfam_matched\trfam_evalue\talignment_score_min\talignment_score_cm\talignment_score_cm_min\t" \
                  "tax_id\tis_bacteria\tgene_id\tgene_location\tgene_description\n"
        if not old_exists:
            all_out_file.write(header)
        for line in all_in_file:
            if skip_lines > 0:
                skip_lines -= 1
                continue
            # 0:design_code, 1:identifier, 2:score, 3:E-value, 4:sequence, 5:round
            # Adding 6:is_rfam_matched	7:rfam_evalue	8:alignment_score_min	9:alignment_score_cm	10:alignment_score_cm_min   11:tax_id 12:is_bacteria  13:gene_id    14:gene_location 15:gene_desctription
            items = line.strip().split('\t')
            design_code = items[0].strip()
            identifier = items[1].strip()
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
            split_ident = compare_rfam.break_id(identifier)
            organism_id = split_ident['code']
            start = split_ident['start']
            end = split_ident['end']
            strand = split_ident['strand']
            subtext = ''
            gene_id = ''
            product = ''
            tax_id = None
            gff_path = get_gff(organism_id)
            if gff_path is not None:
                tax_id = tax_from_gff(gff_path)
                if tax_id is not None:
                    tax_id = result_dive.get_tax_id(organism_id)
                gene_id, product, subtext = nearby_gene(gff_path, int(start), int(end), strand)
            is_bacteria = None
            if tax_id is not None:
                is_bacteria = result_dive.check_ancestor('bacteria', tax_id)
            items.append(str(tax_id))
            items.append(str(is_bacteria))
            items.append(gene_id)
            items.append(subtext)
            items.append(product)
            all_out_file.write('{}\n'.format('\t'.join(items)))
            all_out_file.flush()


def generate_statistics(folder_path: str):
    #generate_common(folder_path)
    extend_final_all(folder_path)


if __name__ == '__main__':
    threshold = [5, 10, 25, 50]
    base_folder = "/opt/home/matan/Dropbox/PHd/RiboSearch/redocm" # final_2018_02_06_11_33_19" #IncManual"
    #for cut in threshold:
    #   left_list = generate_filter(base_folder, cut)
    generate_visuals(base_folder)
    #generate_statistics(base_folder)
