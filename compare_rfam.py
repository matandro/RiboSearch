from typing import List, Dict, Any, TextIO
from Bio import SeqIO
import os
import re


FASTA_LINE_LENGTH = 90
IDENTIFIER_RE = re.compile(r'(?P<code>.+)/(?P<start>\d+)-(?P<end>\d+)(\((?P<strand>[\+|-])\))?')


def break_id(record_id: str) -> Dict[str, Any]:
    matcher = IDENTIFIER_RE.match(record_id)
    if not matcher:
        print("ERROR matching {}".format(record_id))
    code = matcher.group('code')
    start = int(matcher.group('start'))
    end = int(matcher.group('end'))
    strand = matcher.group('strand')
    if strand is None or strand == '':
        strand = '+'
        if end < start:
            strand = '-'
    if strand == '-':
        temp = start
        start = end
        end = temp
    return {"start": start, "end": end, "strand": strand, 'code': code,
            'identifier': "{}/{}-{}({})".format(code,start,end,strand)}


def read_matches(matches_path: str, design_code: str) -> List[Dict[str, Any]]:
    res = []
    with open(matches_path, 'r') as matches_input:
        matches_input.readline()
        for line in matches_input:
            stripped_line = line.strip()
            if stripped_line == '':
                continue
            parts = stripped_line.split('\t')
            if parts[0].strip() != design_code:
                continue
            item = break_id(parts[1])
            item['sequence'] = parts[4]
            res.append(item)
    return res


def read_rfam_fasta(fasta_path: str) -> List[Dict[str, Any]]:
    res = []
    with open(fasta_path, 'rU') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            item = break_id(record.id)
            if item['strand'] == '-':
                item['sequence'] = record.seq.reverse_complement()
            else:
                item['sequence'] = record.seq
            res.append(item)
    return res


def merge_to_rfam(rfams: List[Dict[str, Any]], matches: List[Dict[str, Any]], folder_path: str, prefix:str):
    def test_location(source, target, partial: bool=False):
        if source['strand'] != target['strand']:
            return False
        if source['start'] >= target['start'] and source['end'] <= target['end']:
            return True
        if partial and not source['start'] > target['end'] and not source['end'] < target['start']:
            return True
        return False

    def fasta_add(fasta_writer: TextIO, identifier: str, sequence: str):
        fasta_writer.write("> {}\n".format(identifier))
        for fasta_line in [sequence[i:i + FASTA_LINE_LENGTH] for i in range(0, len(sequence), FASTA_LINE_LENGTH)]:
            fasta_writer.write('{}\n'.format(fasta_line))

    overlaps_path = os.path.join(folder_path, '{}_overlap.fa'.format(prefix))
    match_only_path = os.path.join(folder_path, '{}_match_only.fa'.format(prefix))
    rfam_only_path = os.path.join(folder_path, '{}_rfam_only.fa'.format(prefix))
    old_matches = 0
    rfam_matched = set()
    with open(overlaps_path, 'w') as overlaps_out, open(match_only_path, 'w') as match_only_out,\
        open(rfam_only_path, 'w') as rfam_only_out:
        for match_item in matches:
            is_matched = False
            for rfam_item in rfams:
                if match_item['code'] == rfam_item['code'] and test_location(match_item, rfam_item, partial=True):
                    old_matches += 1
                    rfam_matched.add(rfam_item['identifier'])
                    is_matched = True
                    fasta_add(overlaps_out, match_item['identifier'], match_item['sequence'])
                    break
            if not is_matched:
                fasta_add(match_only_out, match_item['identifier'], match_item['sequence'])
        for rfam_item in rfams:
            if rfam_item['identifier'] not in rfam_matched:
                fasta_add(rfam_only_out, rfam_item['identifier'], rfam_item['sequence'])
    print("Total matched: {}\nTotal Rfam: {}\nOverlap: {}\n".format(len(matches), len(rfams), old_matches))


if __name__ == '__main__':
    BASE_DIR = "/opt/home/matan/Dropbox/PHd/RiboSearch/redocm"
    rfam_records = read_rfam_fasta(os.path.join(BASE_DIR, "RF00167.fa"))
    matches = read_matches(os.path.join(BASE_DIR, 'FINAL_all_ext'), '43_0')
    merge_to_rfam(rfam_records, matches, BASE_DIR, 'cmp')