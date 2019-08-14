from typing import Tuple, List, Dict
from subprocess import Popen
from Bio import SeqIO
import re
import os

FASTA_LINE_LENGTH = 80


def generate_extended(input_list: List[str], out_fasta: str) -> Dict[str, Dict[str, str]]:
    def set_order(start: int, end: int) -> Tuple[int, int]:
        return min(start, end), max(start, end)

    data_map = {}
    expr = re.compile(r'(?P<code>.+)/(?P<start>\d+)-(?P<end>\d+)\((?P<strand>[\+|-])\)')
    for item in input_list:
        matcher = expr.match(item)
        if not matcher:
            print("ERROR matching {}".format(item))
        code = matcher.group('code')
        start, end = set_order(int(matcher.group('start')), int(matcher.group('end')))
        strand = matcher.group('strand')
        data_map[item] = {"start": start, "end": end, "strand": strand, 'code': code}


    for name, data in data_map.items():
        start = data.get('start')
        end = data.get('end')
        code = data.get('code')
        fasta_temp = '{}.fa'.format(code)
        extended_start = max(0, start - 10)
        extended_end = end + 10
        with Popen(['wget', '-O', fasta_temp,
                    'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id={}&from={}&to={}'
                            .format(code, extended_start,extended_end)]) as proc:
            proc.communicate()
            exit_code = proc.wait()
            if exit_code != 0:
                print("ERROR failed to download fasta for {}".format(name))
                continue

        sequence = None
        new_name = None
        try:
            with open(fasta_temp, "rU") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    sequence = record.seq
                    new_name = record.id
                    if strand == '-':
                        sequence = sequence.reverse_complement()
        finally:
            try:
                os.remove(fasta_temp)
            except:
                pass
        data_map[name]['new_name'] = new_name
        data_map[name]['new_sequence'] = sequence

    with open(out_fasta, 'w') as fasta_writer:
        for name, data in data_map.items():
            fasta_writer.write("> {} | {}({})\n".format(name, data['new_name'], data['strand']))
            new_sequence = data['new_sequence']
            for fasta_line in [new_sequence[i:i + FASTA_LINE_LENGTH] for i in range(0, len(new_sequence),
                                                                                    FASTA_LINE_LENGTH)]:
                fasta_writer.write('{}\n'.format(fasta_line))
    return data_map


def read_fasta_ids(fasta_path: str) -> List[str]:
    res = []
    with open(fasta_path, 'rU') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            res.append(record.id)
    return res


def new_old_map(fasta_path: str) -> Dict[str, str]:
    res = {}
    with open(fasta_path, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            res[record.id.split('|')[0].strip()] = record.seq
    return res

if __name__ == "__main__":
    BASE_DIR = "/opt/home/matan/Dropbox/PHd/RiboSearch/redocm"
    item_list = read_fasta_ids(os.path.join(BASE_DIR, 'tranklokase.fa'))
    print(generate_extended(item_list, os.path.join(BASE_DIR, 'extended_tlk.fa')))
