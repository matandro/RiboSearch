import incaRNAtion
import logging
import infernal
from concurrent import futures
from tempfile import NamedTemporaryFile as NTF
import os
import random
import rnafbinv


MAX_WORKER=8
random_gen = random.Random()


def design_seeds(sequence: str, structure: str, amount: int, gc_content: float=0.5):
    return incaRNAtion.run_incaRNAtion(structure, amount, sequence_constraints=sequence, gc_content=gc_content)


def run_design(run_code:int, seed: str, target_sequence: str, target_structure: str, pseudoknots: str=None):
    general_run_logger.info('Starting design {}'.format(run_code))
    result_object = None
    temp_file = None
    try:
        temp_file = NTF(prefix='DESIGN', encoding='utf-8', mode='w', delete=False)
        temp_file.write('TARGET_STRUCTURE={}\n'.format(target_structure))
        temp_file.write('TARGET_SEQUENCE={}\n'.format(target_sequence))
        temp_file.write('TARGET_SEQUENCE={}\n'.format(target_sequence))
        temp_file.write('STARTING_SEQUENCE={}\n'.format(seed))
        temp_file.write('SEED={}\n'.format(random_gen.getrandbits(64)))
        temp_file.write('ITERATION={}\n'.format(1000))
        temp_file.flush()
        temp_file.close()
        result_object = rnafbinv.RNAfbinvCL.main('-f {} --length 5'.format(temp_file.name))
        if result_object is not None:
            if result_object.score < 300 and result_object.score % 100 < 30:
                general_run_logger.info('Finished designing {}'.format(run_code))
                design_logger.info('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(run_code, seed, result_object.sequence,
                                                                       result_object.score, result_object.structure,
                                                                       result_object.bp_dist,
                                                                       result_object.tree_edit_distance))
            else:
                general_run_logger.warning('Finished designing {}, score {}'.format(run_code, result_object.score))
                result_object = None
        else:
            general_run_logger.error('Failed to design sequence. run code: {} seed: {}'.format(run_code, seed))
    except Exception:
        if result_object is None:
            general_run_logger.fatal('Crashed while designing sequences. run code: {} seed: {}'.format(run_code, seed))
    finally:
        if temp_file is not None:
            try:
                os.remove(temp_file.name)
            except:
                pass
    return run_code, result_object


NT_PATH = '/DB/fasta_db/nt/nt'


def run_search(run_code: str, designed_object: rnafbinv.sfb_designer.RnafbinvResult):
    general_run_logger.info('Starting search {}'.format(run_code))
    cm_path = os.path.join(output_dir, '{}.cm'.format(run_code))
    if not infernal.generate_cm(designed_object.sequence, cm_path, designed_object.structure):
        general_run_logger.error('Failed to build covariance model. run code: {}\n{}\n{}'
                                 .format(run_code, designed_object.sequence, designed_object.structure))
        return
    results = infernal.search_cm(cm_path, NT_PATH)
    if results is None:
        general_run_logger.error('Search failed {} {}\n{}'.format(run_code, cm_path, designed_object.sequence))
        return
    general_run_logger.info('Finished search {}, {} results'.format(run_code, len(results)))
    for res_no, res in enumerate(results):
        try:
            sequence = res.get('sequence')
            structure = rnafbinv.vienna.fold(sequence)['MFE']
            res_tree = rnafbinv.shapiro_tree_aligner.get_tree(structure, sequence)
            tree, score = rnafbinv.shapiro_tree_aligner.align_trees(res_tree, target_tree)
            if score < 300 and score % 100 < 30:
                general_run_logger.info('Adding result {}, score {} sequence {}'.format(run_code, score, sequence))
                'seq code\tmatch no\tsequence\tstructure\tscore\ttarget id'
                result_logger.info('{}\t{}\t{}\t{}\t{}\t{}'.format(run_code, res_no, sequence, structure, score,
                                                                   res.get('target name')))
            else:
                general_run_logger.warning('Score too low {} result no {}, score {}, sequence: {}'.format(run_code,
                                                                                                          res_no,
                                                                                                          score,
                                                                                                          sequence))
        except Exception:
            general_run_logger.fatel('Exception in search {}, res no {}, {}'.format(run_code, res_no, res))


def run(sequence: str, structure: str, design_per_seed: int=10, amount_of_seeds: int=500):
    # generate seeds
    seeds = design_seeds(sequence, structure, amount_of_seeds)
    # create executor
    with futures.ThreadPoolExecutor(max_workers=MAX_WORKER) as executor:
        # generate single task for each sequence
        future_sequences = []
        for seed_no, seed in enumerate(seeds):
            for design_no in range(design_per_seed):
                future_sequences.append(executor.submit(run_design, '{}_{}'.format(seed_no, design_no), seed))
        for future in futures.as_completed(future_sequences):
            seq_code, designed_object = future.result()
            if designed_object is not None:
                executor.submit(run_search, seq_code, designed_object)


# Setup new file logger
def setup_logger(logger_name: str, folder: str, formatter: logging.Formatter=logging.Formatter('%(message)s')):
    logger_path = os.path.join(folder, logger_name)
    handler = logging.FileHandler(logger_path)
    handler.setFormatter(formatter)
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    return logger


if __name__ == '__main__':
    # setup vienna location
    rnafbinv.vienna.set_vienna_path('/opt/algorithm/ViennaRNA/bin')
    # setup log files
    output_dir = '/DB/Output/SandD'
    result_logger = setup_logger("match_log", output_dir)
    result_logger.info("seq code\tmatch no\tsequence\tstructure\tscore\ttarget id")
    general_run_logger = setup_logger("error_log", output_dir, logging.Formatter('%(levelname)s::%(asctime)s - %(message)s'))
    design_logger = setup_logger("design_log", output_dir)
    design_logger.info("sec_code\tseed\tsequence\tscore\tstructure\tbp distance\tmotif distance")
    # setup target tree
    PURINE_SEQ = "AGGGUGCCUGAGCCGCGCAUAUAUCGACGGGGAUUCUUCAAAAGCGCCCGCGCGGGGAAACGAAGACGA"
    PURINE_STRUCT = "((((((((...(.(((((.......))))).)........((((((.......))))))..))))))))"
    target_tree = rnafbinv.shapiro_tree_aligner.get_tree(PURINE_STRUCT, PURINE_SEQ)
    # run
    run(PURINE_SEQ, PURINE_STRUCT)