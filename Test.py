from Bio import SeqIO
import sys

test_path = "D:\\Matan\\Dropbox\\PHd\\Michal_Lab\\Dikla\\Analysis\\L_ama_checked_pseudochr.fasta"

with open(test_path, "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        print("record id: {}".format(record.id))
        print("record name: {}".format(record.name))
        print("record description: {}".format(record.description))
        print("record seq: {}".format(record.seq))
        print("record seq: {}".format(record.seq.reverse_complement()))
        print("record seq: {}".format(record.seq.reverse_complement().transcribe()))
        print("record seq[10:20]: {}".format(record.seq[10:20]))
        sys.exit(0)