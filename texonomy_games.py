from typing import Dict, List, Optional

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

###############################################################

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
'''
'''
import os
BASE = "E:"
results = []
sample_dirs = [dir_name for dir_name in os.listdir(BASE) if dir_name.startswith("Sample_") and
               os.path.isdir(os.path.join(BASE, dir_name))]
for sample in sample_dirs:
    single_sample = os.path.join(BASE, sample)
    fastq = [fastq_file for fastq_file in os.listdir(single_sample) if fastq_file.endswith(".fastq.gz")
             and not fastq_file.endswith("clean.fastq.gz")]
    if len(fastq) == 1:
        results.append(os.path.join(single_sample, fastq[0]))
    else:
        print("ERROR: {} : {}".format(sample, fastq))

for item in results:
    print(item)
'''

'''
import result_dive
import os


#BASE_DIR = "/opt/home/matan/Dropbox/PHd/RiboSearch/SandD_finals"
#BASE_DIR = "C:\\Users\\matan\\Dropbox\\PHd\\RiboSearch\\SandD_finals"
BASE_DIR = "D:\\matan\\Dropbox\\PHd\\RiboSearch\\SandD_finals"

HEADER = "0design_code	1identifier	2score	3E-value	4sequence	5is_rfam_matched	6rfam_evalue	" \
         "7alignment_score_min	8alignment_score_cm	9alignment_score_cm_min 10tax_id    11is_bacteria"


with open(os.path.join(BASE_DIR, "FINAL_all_ext_tax")) as input_file:
    input_file.readline()
    for line in input_file:
        words = line.strip().split('\t')
        try:
            tax_id = int(words[10].strip())
            if result_dive.check_ancestor("Viridiplantae", tax_id):
                print(line.strip())
        except ValueError as e:
            print("ERROR: {}".format(e))
        except:
            print("ERROR: {}".format(line.strip()))


matched_rfam = []
matched_location = []
matched_bacteria = []
over_score_rfam = []
total = 0

'''
'''
with open(os.path.join(BASE_DIR, "FINAL_all_ext"), 'r') as input_file, \
    open(os.path.join(BASE_DIR, "FINAL_all_ext_tax"), 'w') as output_file:
    output_file.write("{}\ttax_id\tis_bacteria\n".format(input_file.readline().strip()))
    for line in input_file:
        total += 1
        words = line.strip().split('\t')
        print("Running line {}: {}".format(total, words[1]))
        if words[5] == 1:
            matched_rfam.append(words[0])
            matched_location.append(words[1])
        organism = words[1].split('/')[0]
        org_id = result_dive.get_tax_id(organism)
        is_bacteria = False
        if org_id is not None and result_dive.check_ancestor("bacteria", org_id):
            matched_bacteria.append(words[1])
            is_bacteria = True
        output_file.write("{}\t{}\t{}\n".format(line.strip(), org_id, is_bacteria))
        output_file.flush()

overscored_cm = []
overscored_min = []
THREADHOLD_SCORE = 350
with open(os.path.join(BASE_DIR, "FINAL_all_ext_tax"), 'r') as input_file:
    input_file.readline()
    for line in input_file:
        total += 1
        words = line.strip().split('\t')
        print("Running line {}: {}".format(total, words[1]))
        alignment_score_cm = int(words[8])
        alignment_score_min = int(words[7])
        if alignment_score_cm > THREADHOLD_SCORE:
            overscored_cm.append(words[1])
        if alignment_score_min > THREADHOLD_SCORE:
            overscored_min.append(words[1])
        if words[5] == '1':
            if alignment_score_cm > THREADHOLD_SCORE:
                over_score_rfam.append(words[1])
            matched_rfam.append(words[0])
            matched_location.append(words[1])
        if words[11] == 'True':
            matched_bacteria.append(words[1])

print("Analyzed {} matched, {} overscored cm, {} overscored min".format(total, len(set(overscored_cm)),
                                                                        len(set(overscored_min))))
print("{} matched fit Rfam from {} designed sequences overscores {}".format(len(matched_rfam), len(set(matched_rfam)),
                                                                            len(set(over_score_rfam))))
print("{} matched from bacteria".format(len(matched_bacteria)))

import infernal
BASE_DIR = "/opt/home/matan/Dropbox/PHd/RiboSearch/matches"
'0seq code	1match no	2sequence	3structure	4score	5target id'
total = 0
bacteria = 0
rfam_scan = 0
with open(os.path.join(BASE_DIR, "FINAL_ALL_2018_02_06"), 'r') as input_file,\
    open(os.path.join(BASE_DIR, "FINAL_ALL_2018_02_06_ext"), 'w') as output_file:
    output_file.write("{}\ttax_id\tis_rfam\n".format(input_file.readline().strip()))
    for line in input_file:
        total += 1
        words = line.strip().split('\t')
        target_id = words[1].split('/')[0]
        org_id = result_dive.get_tax_id(target_id)
        seq_db = None
        is_rfam = False
        sequence = words[4]
        try:
            seq_db = infernal.generate_fasta({words[1]: sequence})
            results = infernal.search_cm("/opt/home/matan/Dropbox/PHd/RiboSearch/SandD_finals/RF00167.cm", seq_db.name)
            if results:
                is_rfam = True
                print("Rfamscan: {}".format(line.strip()))
        finally:
            if seq_db is not None and os.path.exists(seq_db.name):
                os.remove(seq_db.name)
        output_file.write("{}\t{}\t{}\n".format(line.strip(), org_id, is_rfam))
'''
'''
import dive_statistics
from collections import OrderedDict
from tempfile import NamedTemporaryFile
BASE_DIR = "/opt/home/matan/Dropbox/PHd/RiboSearch/"
sto_file_name = "RF00167.full.stockholm.txt"

'''
'''
all_file = ""
with open(os.path.join(BASE_DIR, sto_file_name), 'r') as orig_sto,\
    open(os.path.join(BASE_DIR, 'RF00167.full.taxonomy.txt'), 'w') as taxonomy_out:
    taxonomy_out.write('organism\tlocation\ttax_id\n')
    data_map = OrderedDict()
    start_data = False
    max_name = 0
    for line in orig_sto:
        strip_line = line.strip()
        if not start_data:
            if strip_line == "" or line[0] == '#':
                all_file += line
            elif strip_line != "":
                start_data = True

        if start_data and strip_line != "" and strip_line != "//":
            name, data = strip_line.rsplit(maxsplit=1)
            max_name = max(max_name, len(name))
            item = data_map.get(name, "")
            data_map[name] = "{}{}".format(item, data)

    all_file += "\n"
    for name, data in data_map.items():
        all_file += "{} {}\n".format(name.ljust(max_name), data)
        if name[0] != '#':
            organism_loc = name.strip().split('|')[1]
            organism_name, location = organism_loc.strip().split('/')
            print("Running {} {}".format(organism_name, location))
            tax_id = result_dive.get_tax_id(organism_name)
            taxonomy_out.write('{}\t{}\t{}\n'.format(organism_name, location, tax_id))

with NamedTemporaryFile('w', suffix='.sto', dir=BASE_DIR) as ntf:
    ntf.write(all_file)
    ntf.write("//\n")
    ntf.flush()
    dive_statistics.generate_r2r(BASE_DIR, os.path.basename(ntf.name), force=True)



data_list = []
with open(os.path.join(BASE_DIR, 'RF00167.full.taxonomy.txt'), 'r') as taxonomy_in:
    for line in taxonomy_in:
        words = line.strip().split('\t')
        if words[2] == "None":
            tax_id = result_dive.get_tax_id(words[0])
            data_list.append("{}\t{}\t{}\n".format(words[0], words[1], tax_id))
        else:
            data_list.append(line)

with open(os.path.join(BASE_DIR, 'RF00167.full.taxonomy.txt'), 'w') as taxonomy_out:
    for line in data_list:
        taxonomy_out.write(line)
'''

from ete3 import faces, AttrFace, TreeStyle, NodeStyle, NCBITaxa, Tree

ncbi = NCBITaxa()
'''
tax_id_list = []
species_tax_id_map = {}
with open(os.path.join(BASE_DIR, 'RF00167.full.taxonomy.txt'), 'r') as taxonomy_in:
    taxonomy_in.readline()
    for line in taxonomy_in:
        words = line.strip().split('\t')
        tax_id = int(words[2])
        lineage = ncbi.get_lineage(tax_id)
        ranks = ncbi.get_rank(lineage)
        print(ranks)
        species_tax_id = lineage[list(ranks.values()).index('species')]
        tax_id_list.append(tax_id)
        count = species_tax_id_map.get(species_tax_id, 0)
        count += 1
        species_tax_id_map[species_tax_id] = count
        #if not result_dive.check_ancestor('Bacteria', tax_id):
        #    print(line.strip())


def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=30)
        faces.add_face_to_node(N, node, 0, position="aligned")
'''
'''
def get_example_tree():

    # Set dashed blue lines in all leaves
    nst1 = NodeStyle()
    nst1["bgcolor"] = "LightSteelBlue"
    nst2 = NodeStyle()
    nst2["bgcolor"] = "Moccasin"
    nst3 = NodeStyle()
    nst3["bgcolor"] = "DarkSeaGreen"
    nst4 = NodeStyle()
    nst4["bgcolor"] = "Khaki"


    t = Tree("((((a1,a2),a3), ((b1,b2),(b3,b4))), ((c1,c2),c3));")
    for n in t.traverse():
        n.dist = 0

    n1 = t.get_common_ancestor("a1", "a2", "a3")
    n1.set_style(nst1)
    n2 = t.get_common_ancestor("b1", "b2", "b3", "b4")
    n2.set_style(nst2)
    n3 = t.get_common_ancestor("c1", "c2", "c3")
    n3.set_style(nst3)
    n4 = t.get_common_ancestor("b3", "b4")
    n4.set_style(nst4)
    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = False

    ts.mode = "c"
    ts.root_opening_factor = 1
    return t, ts
'''

import result_dive


def count_ancestor(tax_list: List[int], ancestor: str) -> int:
    counter = 0
    for item in tax_list:
        if result_dive.check_ancestor(ancestor, item):
            counter += 1
    return counter


def print_name(tax_list: List[int]):
    x = set()
    for id, name in zip(tax_list, ncbi.translate_to_names(tax_list)):
        if id not in x:
            x.add(id)
            print(name, ':', id)


# ncbi.update_taxonomy_database()
# tree = ncbi.get_topology(species_tax_id_map.keys())
# tree = ncbi.get_topology([13035,13035,13035,13035,1458985,13035,13035,240292,13035,1454205,65093,1085406,1751286,13035,63737,1638788,13035,13035,706587,56107,13035,1638788,13035,2005457,63737,251229,267872,1903187,1973488,1973488,1973488,267872,1903187,706587,1641812,373994,1160,13035,13035,1173020,1454205,1458985])
tree = ncbi.get_topology(
    [8996, 8996, 8996, 8996, 8996, 8996, 8996, 8996, 8996, 8996, 8996, 93934, 93934, 93934, 93934, 93934, 93934, 93934,
     93934, 93934, 93934, 93934, 93934, 8969, 8969, 52644, 52644, 52644, 52644, 9233, 8954, 8954, 8954, 8954, 345164,
     345164, 345164, 345164, 216574, 216574, 216574, 216574, 8839, 8839, 8839, 8839, 8839, 8839, 8839, 8839, 8839, 8839,
     240206, 57412, 9031, 9031, 9031, 9031, 9031, 9031, 9031, 9031, 9031, 9031, 94827, 188379, 121530, 176057, 57068,
     188344, 198806, 198806, 198806, 198806, 198806, 8478, 8478, 9135, 9135, 9135, 9135, 9135, 9135, 9135, 9172, 9172,
     9172, 9172, 9172, 9172, 9172, 13146, 44394, 44394, 44394, 85066, 85066, 85066, 85066, 85066, 85066, 85066, 85066,
     85066, 85066, 85066, 932674, 932674, 932674, 932674, 932674, 932674, 932674, 932674, 932674, 932674, 321398,
     321398, 321398, 321398, 321398, 321398, 321398, 59894, 59894, 59894, 59894, 59894, 59894, 48883, 59729, 13735,
     13735, 13735, 13735, 13735, 181119, 181119, 181119, 181119, 181119, 181119, 9157, 9157, 9157, 9157, 9157, 9157,
     9157, 9157, 38654, 38654, 38654, 38654, 8496, 8496, 8496, 8496, 94835, 8502, 8502, 8502, 8502, 8502, 127582,
     185453, 28737, 28377, 9785, 1230840, 230844, 230844, 230844, 230844, 230844, 230844, 230844, 230844, 230844,
     230844, 230844, 10036, 10116, 10116, 10116, 10116, 10116, 10116, 10116, 10116, 10116, 10116, 10116, 10116, 10116,
     10116, 10116, 10116, 10116, 10116, 10116, 9371, 10020, 9978, 146911, 79684, 79684, 10029, 10029, 10029, 10029,
     10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10093,
     10089, 10089, 10089, 10089, 30608, 10047, 10047, 10047, 10047, 10047, 10047, 10047, 10047, 10047, 10047, 10047,
     10047, 10047, 10047, 10047, 10047, 9986, 9986, 9986, 9986])
# print(tree.get_ascii(show_internal=True, compact=True, attributes=["sci_name"]))
name_tree = Tree(tree.write(features=["sci_name"]))
print(tree.write())

# print_name([1408163, 869754, 1041607])


print_name(
    [8996, 8996, 8996, 8996, 8996, 8996, 8996, 8996, 8996, 8996, 8996, 93934, 93934, 93934, 93934, 93934, 93934, 93934,
     93934, 93934, 93934, 93934, 93934, 8969, 8969, 52644, 52644, 52644, 52644, 9233, 8954, 8954, 8954, 8954, 345164,
     345164, 345164, 345164, 216574, 216574, 216574, 216574, 8839, 8839, 8839, 8839, 8839, 8839, 8839, 8839, 8839, 8839,
     240206, 57412, 9031, 9031, 9031, 9031, 9031, 9031, 9031, 9031, 9031, 9031, 94827, 188379, 121530, 176057, 57068,
     188344, 198806, 198806, 198806, 198806, 198806, 8478, 8478, 9135, 9135, 9135, 9135, 9135, 9135, 9135, 9172, 9172,
     9172, 9172, 9172, 9172, 9172, 13146, 44394, 44394, 44394, 85066, 85066, 85066, 85066, 85066, 85066, 85066, 85066,
     85066, 85066, 85066, 932674, 932674, 932674, 932674, 932674, 932674, 932674, 932674, 932674, 932674, 321398,
     321398, 321398, 321398, 321398, 321398, 321398, 59894, 59894, 59894, 59894, 59894, 59894, 48883, 59729, 13735,
     13735, 13735, 13735, 13735, 181119, 181119, 181119, 181119, 181119, 181119, 9157, 9157, 9157, 9157, 9157, 9157,
     9157, 9157, 38654, 38654, 38654, 38654, 8496, 8496, 8496, 8496, 94835, 8502, 8502, 8502, 8502, 8502, 127582,
     185453, 28737, 28377, 9785, 1230840, 230844, 230844, 230844, 230844, 230844, 230844, 230844, 230844, 230844,
     230844, 230844, 10036, 10116, 10116, 10116, 10116, 10116, 10116, 10116, 10116, 10116, 10116, 10116, 10116, 10116,
     10116, 10116, 10116, 10116, 10116, 10116, 9371, 10020, 9978, 146911, 79684, 79684, 10029, 10029, 10029, 10029,
     10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10093,
     10089, 10089, 10089, 10089, 30608, 10047, 10047, 10047, 10047, 10047, 10047, 10047, 10047, 10047, 10047, 10047,
     10047, 10047, 10047, 10047, 10047, 9986, 9986, 9986, 9986])

items = [8996, 8996, 8996, 8996, 8996, 8996, 8996, 8996, 8996, 8996, 8996, 93934, 93934, 93934, 93934, 93934, 93934,
         93934, 93934, 93934, 93934, 93934, 93934, 8969, 8969, 52644, 52644, 52644, 52644, 9233, 8954, 8954, 8954, 8954,
         345164, 345164, 345164, 345164, 216574, 216574, 216574, 216574, 8839, 8839, 8839, 8839, 8839, 8839, 8839, 8839,
         8839, 8839, 240206, 57412, 9031, 9031, 9031, 9031, 9031, 9031, 9031, 9031, 9031, 9031, 94827, 188379, 121530,
         176057, 57068, 188344, 198806, 198806, 198806, 198806, 198806, 8478, 8478, 9135, 9135, 9135, 9135, 9135, 9135,
         9135, 9172, 9172, 9172, 9172, 9172, 9172, 9172, 13146, 44394, 44394, 44394, 85066, 85066, 85066, 85066, 85066,
         85066, 85066, 85066, 85066, 85066, 85066, 932674, 932674, 932674, 932674, 932674, 932674, 932674, 932674,
         932674, 932674, 321398, 321398, 321398, 321398, 321398, 321398, 321398, 59894, 59894, 59894, 59894, 59894,
         59894, 48883, 59729, 13735, 13735, 13735, 13735, 13735, 181119, 181119, 181119, 181119, 181119, 181119, 9157,
         9157, 9157, 9157, 9157, 9157, 9157, 9157, 38654, 38654, 38654, 38654, 8496, 8496, 8496, 8496, 94835, 8502,
         8502, 8502, 8502, 8502, 127582, 185453, 28737, 28377, 9785, 1230840, 230844, 230844, 230844, 230844, 230844,
         230844, 230844, 230844, 230844, 230844, 230844, 10036, 10116, 10116, 10116, 10116, 10116, 10116, 10116, 10116,
         10116, 10116, 10116, 10116, 10116, 10116, 10116, 10116, 10116, 10116, 10116, 9371, 10020, 9978, 146911, 79684,
         79684, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029, 10029,
         10029, 10029, 10029, 10029, 10093, 10089, 10089, 10089, 10089, 30608, 10047, 10047, 10047, 10047, 10047, 10047,
         10047, 10047, 10047, 10047, 10047, 10047, 10047, 10047, 10047, 10047, 9986, 9986, 9986, 9986]
print(len(items))
print(count_ancestor(items, "Mammalia"))

x = [4932, 559292, 796027, 1344418, 1136231, 294747, 58627, 984485, 573826, 28985, 984485, 983966, 4922, 644223, 981350,
     379508, 871575, 984487, 237561, 237561, 441960, 1229662, 1003335, 4924, 322104, 322104, 4909, 332952, 1170230,
     109264, 441959, 331117, 1509407, 242507, 510516, 330879, 426428, 426428, 426428, 426428, 5061, 425011, 27334,
     1182545, 858893, 653667, 1043004, 5599, 334819, 498257, 1441469, 1286976, 1328760, 91928, 688394, 236234, 1279043,
     1442368, 254056]
for item in x:
    if result_dive.check_ancestor("Aspergillus", item):
        print(item)
