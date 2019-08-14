#!/usr/bin/env python3
'''
Analysis for seed to design results
'''
from Bio import pairwise2
import os
from typing import Tuple, List
import numpy
import concurrent.futures
from rnafbinv import sfb_designer, vienna

# read files
BASE_DIR = "D:\\Matan\\Dropbox\\PHd\\Final Thesis\\data"


def read_single(file_name: str) -> List[Tuple[str, str]]:
    pair_list = []
    with open(os.path.join(BASE_DIR, "{}Seed.txt".format(file_name)), 'r') as input_file:
        input_file.readline()
        for line in input_file:
            words = line.strip().split('\t')
            pair_list.append((words[0].strip('"'), words[1].strip('"')))
    return pair_list


# calc variety
def _calc_single(i: int, j: int, sequence_list: List[str]) -> Tuple[int, int, int]:
    return i, j, pairwise2.align.globalxx(sequence_list[i], sequence_list[j])[0][4]


def calc_socres(sequence_list: List[str]) -> numpy.array:
    scores = numpy.zeros((len(sequence_list), len(sequence_list)))
    futures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        for i in range(0, len(sequence_list) - 1):
            for j in range(i + 1, len(sequence_list)):
                if i != j:
                    futures.append(executor.submit(_calc_single, i, j, sequence_list))
    for future in concurrent.futures.as_completed(futures):
        i, j, score = future.result()
        scores[i][j] = score
        scores[j][i] = score
    return scores


# gather seed and design variety
def calc_varieties(sequence_list_pair: List[Tuple[str, str]]) -> Tuple[numpy.array, numpy.array]:
    seeds = [sequences[0] for sequences in sequence_list_pair]
    designed = [sequences[1] for sequences in sequence_list_pair]
    return calc_socres(seeds), calc_socres(designed)


def _calc_single_diff(i: int, j: int, sequence_list: List[str]) -> Tuple[int, int, int]:
    score = abs(len(sequence_list[i]) - len(sequence_list[j]))
    for index in range(0, min(len(sequence_list[i]), len(sequence_list[j]))):
        if sequence_list[i][index] != sequence_list[j][index]:
            score += 1
    return i, j, score


def calc_diff(sequence_list: List[str]) -> numpy.array:
    scores = numpy.zeros((len(sequence_list), len(sequence_list)))
    futures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        for i in range(0, len(sequence_list) - 1):
            for j in range(i + 1, len(sequence_list)):
                if i != j:
                    futures.append(executor.submit(_calc_single_diff, i, j, sequence_list))
    for future in concurrent.futures.as_completed(futures):
        i, j, score = future.result()
        scores[i][j] = score
        scores[j][i] = score
    return scores


# gather seed and design diffs
def calc_diffs(sequence_list_pair: List[Tuple[str, str]]) -> Tuple[numpy.array, numpy.array]:
    seeds = [sequences[0] for sequences in sequence_list_pair]
    designed = [sequences[1] for sequences in sequence_list_pair]
    return calc_diff(seeds), calc_diff(designed)


sequence_map = {"Random": read_single("rand"), "incaRNAtion": read_single("inca"), "XPT Mutate": read_single("xpt")}
# 1 print sequence varietys
for name, sequence_lists in sequence_map.items():
    print("Calculating diff for {}:".format(name))
    (seed_matrix, designed_matrix) = calc_diffs(sequence_lists)
    numpy.savetxt(os.path.join(BASE_DIR, "{}_diff_seed.txt".format(name)), seed_matrix, delimiter='\t')
    numpy.savetxt(os.path.join(BASE_DIR, "{}_diff_designed.txt".format(name)), designed_matrix, delimiter='\t')
    n = pow(len(sequence_lists), 2)
    print("Variety (seed,design): {},{}".format(sum(sum(seed_matrix)) / n, sum(sum(designed_matrix)) / n))
'''
# 2 print score averages
target_structure = '((((((((...(.(((((.......))))).)........((((((.......))))))..))))))))'
target_sequence = 'NNNNNNNNUNNNNNNNNNNNNNNNNNNNNNNNNUNNNUNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNN'
rnafolder = vienna.LiveRNAfold()
rnafolder.start()
for name, sequence_lists in sequence_map.items():
    print("Calculating scores for {}:".format(name))
    with open(os.path.join(BASE_DIR, "{}_scores.txt".format(name)), 'w') as output_file:
        output_file.write("Seed\tDesigned\tScore\tBP distance\tEdit distance\tEnergy\n")
        for seed, designed in sequence_lists:
            result_object = sfb_designer.RnafbinvResult(designed, {'RNAfold': rnafolder,
                                                                   'target_structure': target_structure,
                                                                   'target_sequence': target_sequence, 'fold': 'MFE'},
                                                        False)
            output_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(seed, designed, result_object.score,
                                                                result_object.bp_dist, result_object.tree_edit_distance,
                                                                result_object.energy))
rnafolder.close()
'''
