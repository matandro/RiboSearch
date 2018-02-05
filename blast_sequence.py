from tempfile import NamedTemporaryFile as NTF
from subprocess import Popen, PIPE
import logging
import os


BLAST_PATH = "/opt/algorithm/blast/bin/"
FASTA_MAX_LENGTH = 80


def generate_query(query_name, sequence):
    query_file = NTF(dir='.', delete=False)
    query_file.write('> {}\n'.format(query_name).encode())
    lines = [sequence[i: i + FASTA_MAX_LENGTH] for i in range(0, len(sequence), FASTA_MAX_LENGTH)]
    for line in lines:
        query_file.write('{}\n'.format(line).encode())
    query_file.close()
    return query_file


def output_blast_analyze(output_file):
    result = []
    with open(output_file, 'r') as blast_res:
        for line in blast_res:
            result.append(line.strip())
    return result


def blast_sequence(sequence, database, query_name='Test', application='blastn', evalue=10, output_str= "6 sseq"):
    results = None
    output_file = None
    sequence_file = None
    try:
        sequence_file = generate_query(query_name, sequence)
        output_file = NTF(dir='.', delete=False)
        param_list = [os.path.join(BLAST_PATH, application), '-db', database,
                      '-query', sequence_file.name, '-evalue', str(evalue), 
                      '-out', output_file.name, '-outfmt', output_str]
        logging.info("starting blast run: {}".format(param_list))
        with Popen(param_list,  stdout=PIPE, stdin=PIPE) as proc:
            proc.communicate()
            if proc.wait() < 0:
                raise Exception(err)
            results = output_blast_analyze(output_file.name)
            logging.info("Finished blast run {} results".format(len(results)))
    except Exception as e:
        logging.error("Failed to run {} on DB {} seq {}. ERROR: {}"
                     .format(application, database, sequence, e))
    finally:
        if sequence_file is not None and os.path.exists(sequence_file.name):
            os.remove(sequence_file.name)
        if output_file is not None and os.path.exists(output_file.name):
            os.remove(output_file.name)
    return results


if __name__ == "__main__":
    sequence = "CAACGGTGCCGTGCGCCTCCTGAGCGTCAAGAACCCGACCATAAGTGAAAAACTTCATAG" \
               "GTGCGGGGAGGGCTGCGCGGCTCTGCGAACCTTCTCGGGCGTCATTTTGTTTTTTATGGC"
    database = "/DB/blast_db/Tbrucei/Tbrucei"
    print(blast_sequence(sequence, database))

