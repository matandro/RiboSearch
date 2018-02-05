import os
import logging
import datetime
import time
import sys
import random
from concurrent import futures

from Bio import Seq

import shapiro_tree_aligner
import incaRNAtion
import rna_designer
import vienna
import infernal
import blast_sequence


MAX_WORKER = 8


def index_single(index_file_path):
    index_map = {}
    with open(index_file_path, 'r') as index_file:
        for line in index_file:
            index, name = line.split(" ",1)
            index_map[name] = index
    return index_map


def index_all(base_folder):
    fasta_index_map = {}
    for db_dir in [a_dir for a_dir in os.listdir(base_folder)
                   if os.path.isdir(os.path.join(base_folder, a_dir))]:
        target_file = os.path.join(base_folder, db_dir, "{}.idx".format(db_dir))
        logging.info("Starting to index {}".format(target_file))
        fasta_index_map[db_dir] = index_single(target_file)
    return fasta_index_map


def gather_fasta_dbs():
    fasta_folder = os.path.join(folder, 'fasta_db')
    all_dbs = [os.path.join(fasta_folder, db_dir) for db_dir in os.listdir(fasta_folder)
               if os.path.isdir(os.path.join(fasta_folder, db_dir))]
    return [os.path.join(one_db, one_db.rsplit('/', 1)[1]) for one_db in all_dbs]


def recover_infernal_sequence(result, fasta_file_path):
    sequence = ''
    db_name = os.path.basename(fasta_file_path)
    key = "{} {}".format(result['target name'], result['description'])
    file_index = indexed_fasta_map.get(db_name)
    if file_index is not None:
        file_index = file_index.get(key)
    if file_index is None:
        # even in error we dont want to lose information
        sequence = '{}_{}_{}({})'.format(key, result['seq from'], result['seq to'], result['strand'])
        logging.warning("Failed to retrieve sequence, key: {}, info: {}".format(key, sequence))
    else:
        with open(fasta_file_path, "r") as fasta_file:
            fasta_file.seek(file_index)
            if result['strand'] == '-':
                start_index = result['seq to']
                end_index = result['seq from']
            else:
                start_index = result['seq from']
                end_index = result['seq to']
            run_to = start_index
            while run_to > 0:
                c = fasta_file.read(1)
                if not c.isspace():
                    run_to -= 1
            run_to = end_index - start_index
            while run_to > 0:
                c = fasta_file.read(1)
                if not c.isspace():
                    sequence += c
                    run_to -= 1
            if result['strand'] == '-':
               seq = Seq(sequence)
               sequence = "{}".format(seq.reverse_complement())
    return sequence


def analyze_res(sequence):
    structure_map = vienna.fold(sequence)
    # generate information from MFE
    structure_mfe = structure_map.get('MFE')
    tree_mfe = shapiro_tree_aligner.get_tree(structure_mfe, sequence)
    match_tree_mfe, match_score_mfe = shapiro_tree_aligner.align_trees(tree_mfe, target_tree)
    # generate information from centroid
    structure_centroid = structure_map.get('centroid')
    tree_centroid = shapiro_tree_aligner.get_tree(structure_centroid, sequence)
    match_tree_centroid, match_score_centroid = shapiro_tree_aligner.align_trees(tree_centroid, target_tree)
    return tree_mfe, match_score_mfe, match_tree_mfe, tree_centroid, match_score_centroid, match_tree_centroid


def add_search_run(source_seq, sequence, seq_code, search_method, db):
    try:
        if sequence != '':
            distance_tuple = analyze_res(sequence)
        else:
            distance_tuple = None, None, None, None, None, None
    except Exception as exc:
        distance_tuple = None, None, None, None, None, None
    result_logger.info("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(seq_code, source_seq, db, sequence,
                                                                             *distance_tuple, search_method))
'''
    "seq code\tdesign sequence\tdb\ttarget sequence\tTree MFE\tdistance MFE\talign tree MFE"
    "\tTree centroid\tdistance centroid\talign tree centroid\tmethod"
'''


def cm_search(sequence, seq_code, fasta_list=None):
    results = []
    # get fasta fils
    if fasta_list is None:
        fasta_dbs = gather_fasta_dbs()
    else:
        fasta_dbs = fasta_list
    # cm build \ calibrate
    cm_path = os.path.join(output_dir, "{}.cm".format(seq_code))
    infernal.generate_cm(sequence, cm_path)
    # search fasta files
    for fasta_file in fasta_dbs:
        single_fasta_res = infernal.search_cm(cm_path, fasta_file)
        for res in single_fasta_res:
            res['file'] = fasta_file 
            #res['sequence'] = recover_infernal_sequence(res, fasta_file)
        results += single_fasta_res
    for res in results:
        add_search_run(sequence, res['sequence'], seq_code, 'cm', res['file'])
    # return results list
    return results


def gather_blast_dbs():
    blast_folder = os.path.join(folder, 'blast_db')
    all_dbs = [os.path.join(blast_folder, db_dir) for db_dir in os.listdir(blast_folder)
               if os.path.isdir(os.path.join(blast_folder, db_dir))]
    return [os.path.join(one_db, one_db.rsplit('/', 1)[1]) for one_db in all_dbs]


def blast_search(sequence, seq_code, db_list=None):
    result = {}
    # run blast seq / dblist
    if db_list is None:
        db_list = gather_blast_dbs()
    for db in db_list:
        result[db] = blast_sequence.blast_sequence(sequence, db)
        logging.info("BLAST search Seq_code: {} found {} results on DB {}".format(seq_code, len(result.get(db, [])), db)) 
    for db, res_list in result.items():
        for res in res_list:
            add_search_run(sequence, res, seq_code, 'blast', db)
    return result


def single_design(run_no, seed, seq_no):
    logging.info("Design - start - run {}, seq_no {}".format(run_no, seq_no))
    # send sequences to RNAfbinv
    designed_sequence = rna_designer.run_rnafbinv(target_structure, target_sequence, iter_no=1000, seed=seed,
                                                  random_seed=random_gen.getrandbits(64))
    designed_structure_map = vienna.fold(designed_sequence)
    logging.info("Design - compare - run {}, seq_no {}".format(run_no, seq_no))
    # generate information from MFE
    designed_structure_mfe = designed_structure_map.get('MFE')
    designed_tree_mfe = shapiro_tree_aligner.get_tree(designed_structure_mfe, designed_sequence)
    match_tree_mfe, match_score_mfe = shapiro_tree_aligner.align_trees(designed_tree_mfe, target_tree)
    # generate information from centroid
    designed_structure_centroid = designed_structure_map.get('centroid')
    designed_tree_centroid = shapiro_tree_aligner.get_tree(designed_structure_centroid, designed_sequence)
    match_tree_centroid, match_score_centroid = shapiro_tree_aligner.align_trees(designed_tree_centroid, target_tree)
    # print info to file
    out_text = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(run_no, seq_no, seed, designed_sequence,
                                                             designed_tree_mfe, match_score_mfe, match_tree_mfe,
                                                             designed_tree_centroid, match_score_centroid,
                                                             match_tree_centroid)
    logging.info("Design - printing results:\n{}".format(out_text))
    design_logger.info(out_text)
    return seq_no, designed_sequence


def single_search(run_no, seq_no, sequence):
    seq_code = "{}_{}".format(run_no, seq_no)
    logging.info("Search - start - seq_code {}".format(seq_code))
    cm_search(sequence, seq_code)
    blast_search(sequence, seq_code)
    logging.info("Search - end - seq_code {}".format(seq_code))
    return True


# Creates a new batch of seeds
def multiple_design(run_no, seeds_per_run, design_per_seed):
    logging.info("Generate seeds run {}".format(run_no))
    # generate seeds
    seeds = incaRNAtion.run_incaRNAtion(target_structure, seeds_per_run, sequence_constraints=target_sequence)
    # create executor
    with futures.ThreadPoolExecutor(max_workers=MAX_WORKER) as executor:
        # generate single task for each sequence
        future_sequences = []
        seq_no = 1
        seed_no = 1
        for seed in seeds:
            for design_no in range(0, design_per_seed):
                future_sequences.append(executor.submit(single_design, run_no, seed, seq_no))
                seq_no += 1
            seed_no +=1
            if seed_no > seeds_per_run:
                break
        for future in futures.as_completed(future_sequences):
            seq_no, designed_sequence = future.result()
            executor.submit(single_search, run_no, seq_no, designed_sequence)

        
def multiple_search(sequence_list):
    with futures.ThreadPoolExecutor(max_workers=MAX_WORKER) as executor:
        future_searches = []
        for sequence in sequence_list:
            future_searches.append(executor.submit(single_search, *sequence))
        for future in futures.as_completed(future_searches):
            try:
                is_good = future.result()
            except Exception as exc:
                logging.error("Error for search: {}".format(exc))


# Setup new file logger
def setup_logger(logger_name, folder, formatter=logging.Formatter('%(message)s')):
    logger_path = os.path.join(folder, logger_name)
    handler = logging.FileHandler(logger_path)
    handler.setFormatter(formatter)
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    return logger


def gather_sequences(design_file_path, match_log_path=None):
    res = None
    searched_sequences = set()
    if match_log_path is not None:
        try:
            with open(match_log_path, 'r') as match_file:
                match_file.readline()
                for line in match_file:
                    searched_sequences.add(line.split('\t', 1)[0].strip())
        except FileNotFoundError:
            pass
    with open(design_file_path, 'r') as design_file:
        res = []
        design_file.readline()
        for line in design_file:
            tokens = [token.strip() for token in line.split('\t')]
            if "{}_{}".format(tokens[0], tokens[1]) not in searched_sequences:
                # collect original run number, original seq number and seqeunce
                res.append((tokens[0], tokens[1], tokens[3]))
    return res


if __name__ == "__main__":
    folder = "/DB/"
    if len(sys.argv) < 4:
        print("Usage: search_runner.py <input sequence> <input structure> <<<amount of runs> <seeds per run> <designed sequence per seed>> | <sequence list>>")
        sys.exit(-1)
    # gather parameters
    target_sequence = sys.argv[1]
    target_structure = sys.argv[2]
    # init random number generator
    random_gen = random.Random()
    # find maximum match
    target_tree = shapiro_tree_aligner.get_tree(target_structure, target_sequence)
    target_score = shapiro_tree_aligner.align_trees(target_tree, target_tree)
    logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(asctime)s:Name[%(name)s]:Thread[%(thread)d] - %(message)s')
    output_dir = os.path.join(folder, "Output")
    # new method infernal has it's out index
    #indexed_fasta_map = index_all(os.path.join(folder, "fasta_db"))
    # Start runs
    if len(sys.argv) == 6:
        # init log (result) files
        run_time_stamp = datetime.datetime.fromtimestamp(time.time()).strftime("%Y_%m_%d_%H_%M_%S")
        output_dir = os.path.join(output_dir, run_time_stamp)
        os.mkdir(output_dir)
        result_logger = setup_logger("match_log", output_dir)
        result_logger.info("seq code\tdesign sequencei\tdb\ttarget sequence\tTree MFE\tdistance MFE\talign tree MFE"
                       "\tTree centroid\tdistance centroid\talign tree centroid\tmethod")
        design_logger = setup_logger("design_log", output_dir)
        design_logger.info("run\tseq no\tseed\tsequence\tTree MFE\tdistance MFE\talign tree MFE"
                           "\tTree centroid\tdistance centroid\talign tree centroid")
        for i in range(1, int(sys.argv[3]) + 1):
            multiple_design(i, int(sys.argv[4]), int(sys.argv[5]))
    else:
        design_logger = logging.getLogger('dummy') # dummy logger for design since we might no need one
        output_dir = os.path.dirname(sys.argv[3])
        match_log_path = os.path.join(output_dir, "match_log")
        if not os.path.exists(match_log_path):
            match_log_path = None
        sequence_list = gather_sequences(sys.argv[3], match_log_path)
        result_logger = setup_logger("match_log", output_dir) # appending to exising file
        if match_log_path is None:
            result_logger.info("seq code\tdesign sequencei\tdb\ttarget sequence\tTree MFE\tdistance MFE\talign tree MFE"
                               "\tTree centroid\tdistance centroid\talign tree centroid\tmethod")
        multiple_search(sequence_list)

