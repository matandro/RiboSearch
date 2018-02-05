import os
import logging
from tempfile import NamedTemporaryFile as NTF
import itertools
import sys
from blast_sequence import blast_sequence

import infernal
import shapiro_tree_aligner


def is_novel(sequence, cm_file):
    res = True
    fasta_file = None
    try:
        fasta_file = NTF(dir='.', delete=False)
        fasta_file.write('> seq\n')
        for item in map(''.join, itertools.zip_longest(*[iter(sequence)]* 80, fillvalue='')):
            fasta_file.write('{}\n'.format(item))
        fasta_file.close()
        cm_res = infernal.search_cm(cm_file, fasta_file.name)
        if cm_res is not None and len(cm_res) > 0:
            res = False
    finally:
        if fasta_file is not None and os.path.exists(fasta_file.name):
            os.remove(fasta_file.name)
    return res


def remove_existing(match_folder, cm_file, filter_score):
    def get_index(headers, header_name, def_value):
        try:
            index = headers.index(header_name)
        except ValueError:
            logging.warning('Could not find "{}" in header. assuming position {}'.format(header_name, def_value))
            index = def_value
        return index
    def gather_blast(db):
        short_db = db.rsplit('/', 1)
        return "/DB/blast_db/{}/{}".format(short_db, short_db)
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
            score_index = get_index(header, "distance centroid", 8)
            db_index = get_index(header, "db", 2)
            out_file.write("{}\n".format(header))
            for line in in_file:
                if line.strip() == '' or line[:1] == '#':
                    continue
                values = [value.strip() for value in line.strip().split('\t')]
                sequence = values[sequence_index]
                score = float(values[score_index])
                database = gather_blast(values[db_index])
                if score < filter_score and is_novel(sequence, cm_file):
                    sequence_info = blast_sequence(sequence, database, output_str="6 pident sacc sstart send")
                    # need to decide how to analyze results and what should be added
                    collapsed_info = test_info(sequence_info)
                    out_file.write("{}\n".format(line))


def analyze_all(output_folder, cm_file, filter_score):
    subdirs = [outdir for outdir in os.listdir(output_folder)
               if os.path.isdir(outdir)]
    for outdir in subdirs:
        remove_existing(outdir, cm_file)


def get_max_score(seq, struct):
    tree = shapiro_tree_aligner.get_tree(struct, seq)
    _, score = shapiro_tree_aligner.align_trees(tree, tree)
    return score


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    if len(sys.argv) > 5:
        folder = sys.argv[1]
        cm_file = sys.argv[2]
        sequence = sys.argv[3]
        structure = sys.argv[4]
        min_score = float(sys.argv[5])
    else:
        logging.warning("""
                        No arguments recieved: filter_results.py <output folder> <cm file> <sequence> <structure> <min score%>
                        Example: python3 filter_results.py '/DB/Output/' 'purine.cm' 'NNNNNNNNUNNNNNNNNNNNNNNNNNNNNNNNNUNNNUNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNN' '((((((((...(.(((((.......))))).)........((((((.......))))))..))))))))' 60.0
                        Using defaults
                        """)
        folder = '/DB/Output/'
        cm_file = 'purine.cm'
        sequence = 'NNNNNNNNUNNNNNNNNNNNNNNNNNNNNNNNNUNNNUNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNN'
        structure = '((((((((...(.(((((.......))))).)........((((((.......))))))..))))))))'
        min_score =  60.0
    max_score = get_max_score(sequence, structure)
    min_score = max_score * (min_score / 100.0)
    analyze_all(folder, cm_file, min_score)

