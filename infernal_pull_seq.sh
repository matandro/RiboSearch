#!/bin/sh

# Usage: sh infernal_pull_seq.sh <tlbout file> <input fasta file> <output fasta>

TLBOUT_FILE=$1
INPUT_FASTA=$2
OUTPUT_FASTA=$3

ESL_SFETCH_LOC="/opt/algorithm/infernal/bin/esl-sfetch"

if [ -s $TLBOUT_FILE ]; then
    cat $TLBOUT_FILE | grep -v ^\# | awk '{ printf("%s/%d-%d %d %d %s\n", $1, $8, $9, $8, $9, $1); }' | $ESL_SFETCH_LOC \-Cf $INPUT_FASTA - > $OUTPUT_FASTA
fi