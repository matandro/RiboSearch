from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from tempfile import NamedTemporaryFile as ntf
import os
import logging
import sys


FASTA_MAX_LENGTH = 80
BLAST_PROG = 'blastn'


def create_ss_fasta(sequence, topic='seq'):
    fasta_file = ntf(prefix='blast', delete=False, mode='w')
    fasta_file.write('> {}\n'.format(topic).encode())
    lines = [sequence[i: i + FASTA_MAX_LENGTH] for i in range(0, len(sequence), FASTA_MAX_LENGTH)]
    for line in lines:
        fasta_file.write('{}\n'.format(line).encode())
        fasta_file.close()
    return fasta_file.name


def create_ms_fasta(seq_topic_map):
    fasta_file = ntf(prefix='blast', delete=False, mode='w')
    for sequence, topic in seq_topic_map.items():
        fasta_file.write('> {}\n'.format(topic).encode())
        lines = [sequence[i: i + FASTA_MAX_LENGTH] for i in range(0, len(sequence), FASTA_MAX_LENGTH)]
        for line in lines:
            fasta_file.write('{}\n'.format(line).encode())
    fasta_file.close()
    return fasta_file.name


def basic_parser(blast_record):
    res = {'name': blast_record.Header.query,
           'db': blast_record.Header.database}
    alignment_list = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            single_match = {'title': alignment.title,
                            'hit_id': alignment.hit_id,
                            'score': hsp.score,
                            'bits': hsp.bits,
                            'expect': hsp.expect,
                            'length': hsp.align_length,
                            'identities': hsp.identities,
                            'positives': hsp.positives,
                            'strand': hsp.strand,
                            'query': hsp.query,
                            'query_start': hsp.query_start,
                            'query_end': hsp.query_end,
                            'match': hsp.match,
                            'sbjct': hsp.sbjct,
                            'sbjct_start': hsp.sbjct_start,
                            'sbjct_end': hsp.sbjct_end}
            alignment_list.append(single_match)
    res['alignments'] = alignment_list
    return res


def full_alignment_parser(blast_record):
    res = basic_parser(blast_record)
    alignment_list = []
    for alignment in res['alignments']:
        if alignment['identities'] == alignment['length']:
            alignment_list.append(alignment)
    res['alignments'] = alignment_list
    return res


def blast_file(fasta_path, blast_db='nt', parser=basic_parser):
    results = []
    record = SeqIO.read(fasta_path).read()
    result_handle = NCBIWWW.qblast(BLAST_PROG, blast_db, record.seq, megablast=True)
    blast_records = NCBIXML.parse(result_handle)
    for single_record in blast_records:
        # each run is a single sequence search from fasta_path
        results.append(parser(single_record))
    return results


MAX_GROUP_SIZE = 50


def run_filtered_file(file_path, out_path):
    group = []
    with open(file_path, 'r') as in_file, open(out_path, 'w') as out_file:
        out_file.write('seq_id\tsequence\thit id\tstrand\tstart\tend\ndescription\t')
        in_file.readline()
        for line in in_file:
            if line.strip() == '' or line is None:
                continue
            values = line.split('\t')
            search_id = values[0]
            sequence = values[3]
            group.append((search_id, sequence))
            if len(group) == MAX_GROUP_SIZE:
                run_group(group, out_file)
                group = []
        if len(group) > 0:
            run_group(group, out_file)


def run_group(sequence_list, out_file):
    fasta_temp = None
    try:
        sequence_map = {}
        for seq_index in range(0, len(sequence_list)):
            sequence_map['SEQ{}'.format(seq_index)] = sequence_list[seq_index][2]
        fasta_temp = create_ms_fasta(sequence_map)
        results = blast_file(fasta_temp, parser=full_alignment_parser)
        ret_val = [[]] * len(sequence_list)
        for res in results:
            try:
                seq_index = int(res['name'].split('SEQ')[1].strip())
            except ValueError:
                logging.error("Could not analyze sequence index res map:{}".format(res))
            search_id, sequence = sequence_list[seq_index]
            for alignment in res.alignments:
                out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(search_id, sequence,
                                                                     alignment['hit_id'],
                                                                     alignment['strand'],
                                                                     alignment['sbjct_start'],
                                                                     alignment['sbjct_end'],
                                                                     alignment['title'],))
    finally:
        if fasta_temp is not None:
            os.remove(fasta_temp)


if __name__ == "__main__":
    run_filtered_file(sys.argv[1], sys.argv[2])