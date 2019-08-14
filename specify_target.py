from rnafbinv import vienna, shapiro_tree_aligner
from tempfile import NamedTemporaryFile
from typing import Tuple, Dict, List, Optional
import dive_statistics
import increased_pull
import result_dive
import infernal
import shutil
import os


def read_sto(sto_file: str) -> Dict[str, Tuple[str, str]]:
    structure = ''
    seq_map = {}
    with open(sto_file, 'r') as sto_input:
        for line in sto_input:
            stripped_line = line.strip()
            if stripped_line == "" or stripped_line == '//':
                continue
            if line.startswith("#=GC SS_cons"):
                structure += stripped_line.rsplit(maxsplit=1)[1]
            elif line[0] != "#":
                name, sequence = stripped_line.split(maxsplit=1)
                seq_map[name] = seq_map.get(name, "") + sequence
    for name, sequence in seq_map.items():
        final_seq = ''
        final_struct = ''
        if len(sequence) != len(structure):
            print("ERROR, BAD LENGTH\n{}\n{}".format(sequence, structure))
            exit(-3)
        origin_map_for = {}
        origin_map_rev = {}
        origin_new_map = {}
        new_origin_map = {}
        bracket_stack = []
        index = 0
        for tseq_char, struct_char in zip(sequence, structure):
            if struct_char in '[{(<':
                bracket_stack.append(index)
            elif struct_char in ']})>':
                close_index = bracket_stack.pop()
                origin_map_for[close_index] = index
                origin_map_rev[index] = close_index
            if tseq_char.upper() in 'AGCU':
                if struct_char in '[]{}()<>':
                    origin_new_map[len(final_struct)] = index
                    new_origin_map[index] = len(final_struct)
                final_struct += struct_char
                final_seq += tseq_char.upper()
            index += 1
        fix_structure = ''
        for index, struct_char in enumerate(final_struct):
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
        final_struct = fix_structure
        seq_map[name] = (final_seq, final_struct)
    return seq_map


def clean_sto(removed_names: str, sto_name: str):
    with NamedTemporaryFile('w') as named_input:
        with open(sto_name, 'r') as sto_input:
            for line in sto_input:
                remove = False
                for name in removed_names:
                    if name in line:
                        remove = True
                        break
                if not remove:
                    named_input.write(line)
        named_input.flush()
        shutil.copy(named_input.name, sto_name)


def recreate_cm(folder_path: str, design_code: str, desc_keywords: List[str], test_name: str, sequence: str,
                structure: str, identifier_keywords: List[str]= None, new_seq_map: Dict[str, str] = {},
                top_taxonomy: str=None, score_cutoff: int=250):
    def keyword_relevant(check_item: str, keymap: List[str]) -> bool:
        if keymap is None:
            return True
        for key in keymap:
            if key.lower() in check_item.lower():
                return True
        return False

    def check_taxonomy(tax_id: Optional[int]) -> bool:
        if top_taxonomy is None:
            return True
        elif tax_id is None:
            return False
        return result_dive.check_ancestor(top_taxonomy, tax_id)

    def check_c_or_u(named_tree) -> Tuple[bool, bool]:
        _, score = shapiro_tree_aligner.align_trees(named_tree, c_tree)
        check_c = score < 1000
        _, score = shapiro_tree_aligner.align_trees(named_tree, u_tree)
        check_u = score < 1000
        return check_c, check_u

    tree = shapiro_tree_aligner.get_tree(
        ".....((((((((...(.(((((.......))))).)........((((((.......))))))..)))))))).....",
        "NNNNNNNNNNNNNUNNNNNNNNNNNNNNNNNNNNNNNNUNNNUNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNNNNNNN")
    c_tree = shapiro_tree_aligner.get_tree(
        ".....((((((((...(.(((((.......))))).)........((((((.......))))))..)))))))).....",
        "NNNNNNNNNNNNNUNNNNNNNNNNNNNNNNNNNNNNNNUNNNUNNNNNNNNNNNNNNNNNNNNNNCNNNNNNNNNNNNN")
    u_tree = shapiro_tree_aligner.get_tree(
        ".....((((((((...(.(((((.......))))).)........((((((.......))))))..)))))))).....",
        "NNNNNNNNNNNNNUNNNNNNNNNNNNNNNNNNNNNNNNUNNNUNNNNNNNNNNNNNNNNNNNNNNUNNNNNNNNNNNNN")
    gather_results = {}
    full_results = {}
    with open(os.path.join(folder_path, "FINAL_all_ext"), 'r') as input_all:
        input_all.readline()
        for line in input_all:
            if line.strip() == '':
                continue
            parts = line.strip().split('\t')
            try:
                tax_id = int(parts[11])
            except ValueError:
                tax_id = None
            if parts[0].strip() == design_code and keyword_relevant(parts[15], desc_keywords) and\
                    keyword_relevant(parts[1], identifier_keywords) and check_taxonomy(tax_id):
                read_sequence = new_seq_map.get(parts[1], parts[4])
                gather_results[parts[1]] = read_sequence
                full_results[parts[1]] = {'sequence': read_sequence, 'rfam_eval': parts[7], 'tax_id': parts[11],
                                          'gene_id': parts[13], 'gene_loc': parts[14], 'gene_desc': parts[15]}
    if len(gather_results) == 0:
        print("No results from design {} match keywords {}".format(design_code, desc_keywords))
        exit(0)
    try:
        cm = NamedTemporaryFile('w', suffix=".cm")
        cm_file_name = cm.name
        cm.close()
        sto = NamedTemporaryFile('w', suffix=".sto")
        sto_file_name = sto.name
        sto.close()
        real_cm = os.path.join(folder_path, "{}.cm".format(test_name))
        if not os.path.exists(real_cm):
            if not infernal.generate_single_seq_cm(sequence, cm_file_name, structure):
                print("Failed to generate single sequence-structure cm file for {}".format(design_code))
                exit(-1)
            shutil.copy(cm_file_name, real_cm)
        if not infernal.align_sequences(gather_results, os.path.join(folder_path, real_cm), sto_file_name):
            print("Failed to generate alignment for {}".format(design_code))
            exit(-2)
        removed_names = set()
        cm_map = read_sto(sto.name)
        with open(os.path.join(folder_path, "{}.txt".format(test_name)), 'w') as tbl_out:
            tbl_out.write('name\tscore\ttax id\trfam evalue\tgene id\tgene location\tgene description\tcm sequence\t'
                          'cm structure\toriginal sequence\tguanine_bind\tadenine_bind\n')
            for name, (cmsequence, cmstructure) in cm_map.items():
                named_tree = shapiro_tree_aligner.get_tree(cmstructure, cmsequence)
                _, score = shapiro_tree_aligner.align_trees(named_tree, tree)
                is_c, is_u = check_c_or_u(named_tree)
                if score > score_cutoff:
                    removed_names.add(name)
                    print("Removing {} score {}\n{}\n{}".format(name, score, cmsequence, cmstructure))
                else:
                    res_item = full_results.get(name)
                    tbl_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(name, score, res_item.get('tax_id'),
                                                                                    res_item.get('rfam_eval'),
                                                                                    res_item.get('gene_id'),
                                                                                    res_item.get('gene_loc'),
                                                                                    res_item.get('gene_desc'),
                                                                                    cmsequence, cmstructure,
                                                                                    res_item.get('sequence'),
                                                                                    is_c, is_u))
        clean_sto(removed_names, sto_file_name)
        real_sto = os.path.join(folder_path, "{}.sto".format(test_name))
        shutil.copy(sto_file_name, real_sto)
        dive_statistics.generate_r2r(folder_path, real_sto, force=True)

    finally:
        try:
            os.remove(cm_file_name)
        except:
            pass
        try:
            os.remove(sto_file_name)
        except:
            pass



if __name__ == "__main__":
    BASE_DIR = "/opt/home/matan/Dropbox/PHd/RiboSearch/redocm"

    # Transketolase
    sequence = "TCTATCGGGCTTGGTGAGGATGGTCCAACGCACCAACCTATTGAAACTCTGGCTCACTTGAGGGCTATTCCAACCCGGCATGT"
    structure = "....(((((.((((..(((.((((.......)))).)))..........((((((......))))))..)))))))))....."
    DESIGN = "M_1"
    KEYWORDS = ['transketolase']#'tkl',
    top_taxonomy = "fungi"
    #top_taxonomy = None
    name = "transketolase_ext_fungi"
    extended_sequences = increased_pull.new_old_map(os.path.join(BASE_DIR, 'extended_tlk.fa'))
    '''
    # mammalian
    sequence = "AGACGAGAUUUUAACAUCUGCAAAUUGAGGUAGUUUUAUUUCUUCCUUUCCAUUUUGGAAGCCUUUUGUCU"
    structure = "((((((((....(((((((........)))).)))......(((((..........)))))..))))))))"
    DESIGN = "M_6"
    KEYWORDS = None
    extended_sequences = {} #increased_pull.new_old_map(os.path.join(BASE_DIR, 'M_6.fa'))
    top_taxonomy = "mammalian"
    name = "M_6_mammalian"
    

    sequence = "CUCCCGUUACCGUCGCUGGUUAAGCGCAGGUUUGUAGAACUCGGAUAUAAAAAUCCAUCCGGGGG"
    structure = "((((((..(((..((((.....))))..)))...........((((......))))...))))))"
    DESIGN = "156_2"
    KEYWORDS = ['myosin']
    extended_sequences = {} #increased_pull.new_old_map(os.path.join(BASE_DIR, 'M_6.fa'))
    top_taxonomy = 'mammalia'
    name = "156_2_myosin_mammalia"
'''
    recreate_cm(BASE_DIR, DESIGN, KEYWORDS,name, sequence, structure, new_seq_map=extended_sequences,
                top_taxonomy=top_taxonomy, score_cutoff=250)
