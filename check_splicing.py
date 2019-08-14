from subprocess import Popen
from Bio import SeqIO
import compare_rfam
import tempfile
import gffutils
import os
import result_dive


def get_sequence(gff_folder: str, name: str, gene_id: str):
    sequence = None
    db_file = os.path.join(gff_folder, "{}.db".format(name))
    gffdb = gffutils.FeatureDB(db_file)
    item = gffdb[gene_id]
    start = item.start
    end = item.end
    strand = item.strand
    fasta_temp = tempfile.NamedTemporaryFile('w', suffix=".fa", delete=False)
    fasta_temp.close()
    read = False
    try:
        with Popen(['wget', '-O', fasta_temp.name,
                    'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id={}&from={}&to={}'
                            .format(name, start, end)]) as proc:
            proc.communicate()
            exit_code = proc.wait()
            if exit_code != 0:
                print("ERROR failed to download fasta for {}".format(name))
            else:
                read = True
        if read:
            with open(fasta_temp.name, "rU") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    sequence = record.seq
                    if strand == '-':
                        sequence = sequence.reverse_complement()
    finally:
        try:
            os.remove(fasta_temp.name)
        except:
            pass
    return sequence


if __name__ == "__main__":
    BASE_DIR = "/opt/home/matan/Dropbox/PHd/RiboSearch"
    gff_dir = os.path.join(BASE_DIR, 'gff')
    tsv_dir = os.path.join(BASE_DIR, 'redocm')
    TSV_FILE = "156_2_myosin.txt" #"extended_tlk_fungi.txt"
    with open(os.path.join(tsv_dir, TSV_FILE), 'r') as input_tsv,\
        open(os.path.join(tsv_dir, "{}.fa".format(TSV_FILE[:-4])), 'w') as fasta_out:
        input_tsv.readline()
        for line in input_tsv:
            parts = line.strip().split('\t')
            if len(parts) == 1: # or not result_dive.check_ancestor('Homo sapiens',int(parts[2])):
                continue
            item_broken = compare_rfam.break_id(parts[0])
            gene_id = parts[4][1:-1].strip("'")
            sequence = get_sequence(gff_dir, item_broken['code'], gene_id)
            location = sequence.transcribe().find(parts[7])
            fasta_out.write("> {} as {}-{}\n".format(item_broken['identifier'], location, location+len(parts[7])))
            for fasta_line in [sequence[i:i + 80] for i in range(0, len(sequence), 80)]:
                fasta_out.write('{}\n'.format(fasta_line))
