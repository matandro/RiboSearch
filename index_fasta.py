import sys
import os
import logging


def index_single(db_fasta_file):
    with open(db_fasta_file, 'r') as fasta_file:
        with open("{}.idx".format(db_fasta_file), 'w') as out_file:
            is_name = False
            name = ""
            while True:
                c = fasta_file.read(1)
                if not c:
                    break # EOF
                if c == ">" and not is_name:
                    is_name = True
                elif is_name:
                    if c == "\n":
                        out_file.write("{}\t{}\n".format(fasta_file.tell(), name))
                        name = ""
                        is_name = False
                    else:
                        name += c


def index_all(base_folder):
    for db_dir in [a_dir for a_dir in os.listdir(base_folder) 
                   if os.path.isdir(os.path.join(base_folder, a_dir))]:
        target_file = os.path.join(base_folder, db_dir, db_dir)
        logging.info("Starting to index {}".format(target_file))
        index_single(target_file)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    folder = "/DB/fasta_db/"
    if len(sys.argv) > 1: 
        folder = sys.argv[1]
    index_all(folder)

