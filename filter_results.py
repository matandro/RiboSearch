import os
import logging
from tempfile import NamedTemporaryFile as NTF
import itertools
import sys
from blast_sequence import blast_sequence

import infernal
import shapiro_tree_aligner
import search_runner


def is_novel(sequence, cm_file):
    res = True
    fasta_file = None
    try:
        fasta_file = NTF(dir='.', mode='w', delete=False)
        fasta_file.write('> seq\n')
        for item in map(''.join, itertools.zip_longest(*[iter(sequence)]* 80, fillvalue='')):
            fasta_file.write('{}\n'.format(item))
        fasta_file.close()
        cm_res = infernal.search_cm(cm_file, fasta_file.name, res_type=infernal.ResType.TLBOUT)
        if cm_res is not None and len(cm_res) > 0:
            res = False
    finally:
        if fasta_file is not None and os.path.exists(fasta_file.name):
            os.remove(fasta_file.name)
    return res


def search_blast():
    def gather_blast(db):
        short_db = db.rsplit('/', 1)[1]
        return "/DB/blast_db/{}/{}".format(short_db, short_db)
    database = gather_blast(values[db_index])
    sequence_info = blast_sequence(sequence, database, output_str="6 sacc sstart send sstrand", extra_options=['-perc_identity', '100'])
    # need to decide how to analyze results and what should be added
    collapsed_info = test_info(sequence_info)

 
'''    
seq code\tdesign sequence\tdb\ttarget sequence\tTree MFE\tdistance MFE\talign tree MFE
                               \tTree centroid\tdistance centroid\talign tree centroid\tmethod 
'''
def remove_existing(match_folder, cm_file, filter_score, recalc_res=False):
    def get_index(headers, header_name, def_value):
        try:
            index = headers.index(header_name)
        except ValueError:
            logging.warning('Could not find "{}" in header. assuming position {}'.format(header_name, def_value))
            index = def_value
        return index
    def test_info(info_list):
        real_list = [",".join(item.split('\t')) for item in info_list]
        return ";".join(real_list)
    match_file_path = os.path.join(match_folder, "match_log")
    if not os.path.exists(match_file_path):
        logging.warning("No match_log file in {}".format(match_folder))
    else:
        with open("{}_nocm".format(match_file_path), 'w') as out_file, open(match_file_path, 'r') as in_file:
            header = [item.strip() for item in in_file.readline().strip().split('\t')]
            sequence_index = get_index(header, "target sequence", 3)
            db_index = get_index(header, "db", 2)
            out_file.write("{}\n".format(header))
            for line in in_file:
                if line.strip() == '' or line[:1] == '#':
                    continue
                values = [value.strip() for value in line.strip().split('\t')]
                sequence = values[sequence_index]
                logging.info("Running sequence: {}".format(sequence))
                if recalc_res:
                    tree_mfe, match_score_mfe, match_tree_mfe, tree_centroid, match_score_centroid, match_tree_centroid = search_runner.analyze_res(sequence, target_tree)
                    score = match_score_centroid
                    values[get_index(header, "Tree MFE", 4)] = tree_mfe
                    values[get_index(header, "distance MFE", 5)] = match_score_mfe
                    values[get_index(header, "align tree MFE", 6)] = match_tree_mfe
                    values[get_index(header, "Tree centroid", 7)] = tree_centroid
                    values[get_index(header, "distance centroid", 8)] = match_score_centroid 
                    values[get_index(header, "align tree centroid", 9)] = match_tree_centroid
                    line = "{}\n".format("\t".join([str(item) for item in values]))
                else:
                    score = float(values[get_index(header, "distance centroid", 8)])
                if score > filter_score and is_novel(sequence, cm_file):
                    logging.info("LINE ADDED: [{}\n], score: {}".format(line, score))
                    out_file.write("{}\n".format(line))


def analyze_all(output_folder, cm_file, filter_score, recalc_res):
    subdirs = [os.path.join(output_folder, outdir) for outdir in os.listdir(output_folder)
               if os.path.isdir(os.path.join(output_folder, outdir))]
    for outdir in subdirs:
        logging.info("Running single analysis for {}".format(outdir))
        remove_existing(outdir, cm_file, filter_score, recalc_res)


def get_max_score(seq, struct):
    tree = shapiro_tree_aligner.get_tree(struct, seq)
    _, score = shapiro_tree_aligner.align_trees(tree, tree)
    return tree, score


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    if len(sys.argv) > 5:
        folder = sys.argv[1]
        cm_file = sys.argv[2]
        target_sequence = sys.argv[3]
        target_structure = sys.argv[4]
        min_score = float(sys.argv[5])
        recalc = False
        if len(sys.argv) > 6:
            recalc = (sys.argv[6].strip().lower() == 'true')
    else:
        logging.warning("""No arguments recieved: filter_results.py <output folder> <cm file> <sequence> <structure> <min score%> [recalc]
Example: python3 filter_results.py '/DB/Output/' 'purine.cm' 'NNNNNNNNUNNNNNNNNNNNNNNNNNNNNNNNNUNNNUNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNN' '((((((((...(.(((((.......))))).)........((((((.......))))))..))))))))' 85.0 True
Using defaults""")
        folder = '/DB/Output/'
        cm_file = 'purine.cm'
        target_sequence = 'NNNNNNNNUNNNNNNNNNNNNNNNNNNNNNNNNUNNNUNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNN'
        target_structure = '((((((((...(.(((((.......))))).)........((((((.......))))))..))))))))'
        min_score =  95.0
        recalc= True
    target_tree, max_score = get_max_score(target_sequence, target_structure)
    logging.info("Optimal score is: {}".format(max_score))
    min_score = max_score * (min_score / 100.0)
    analyze_all(folder, cm_file, min_score, recalc)

