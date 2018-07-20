'''
import rpy2
print(rpy2.__version__)

from rpy2.rinterface import R_VERSION_BUILD
print(R_VERSION_BUILD)


import rpy2.robjects as robjects

BAYESIAN_CLUSTER_DIR = 'D:/Matan/Dropbox/PHd/Michal_Lab/Limor DBscan/baysian_cluster/'
BAYESIAN_CLUSTER_PATH = 'Bayesian_analysis_alpha.R'


r_setwd = robjects.r['setwd']
print(r_setwd(BAYESIAN_CLUSTER_DIR))

r_getwd = robjects.r['getwd']
print(r_getwd())


r_source = robjects.r['source']
print(r_source(BAYESIAN_CLUSTER_PATH))


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
'''

import re


def extract_gene(mstr):
    res = None
    match = re.match(r'.*gene=(?P<gene_id>[^\s]+).*', mstr)
    if match is not None:
        res = match.group('gene_id')
    else:
        print('failed to match line: {}'.format(mstr))
    return res

infile = 'D:\\Matan\\Dropbox\\PHd\\Michal_Lab\\Rohit Analysis\\ms\\MaxQuant\\fasta_header.txt'
outfile = 'D:\\Matan\\Dropbox\\PHd\\Michal_Lab\\Rohit Analysis\\ms\\MaxQuant\\fasta_gene.txt'

header_set = set()
with open(outfile, 'w') as outf:
    with open(infile, 'r') as inf:
        for line in inf:
            fasta_line = line.strip('\n')
            if fasta_line not in header_set:
                header_set.add(fasta_line)
                gene_name = extract_gene(fasta_line)
                outf.write('{}\t{}\n'.format(fasta_line, '' if gene_name is None else gene_name))

