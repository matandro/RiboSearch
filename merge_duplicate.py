import os


def rem_dups(in_file_path, out_file_path):
    with open(in_file_path, 'r') as in_file, open(out_file_path, 'w') as out_file:
        # write header
        out_file.write("{}".format(in_file.readline()))
        # gather and merge lines
        line_map = {}
        for line in in_file:
            seq_code, rest_of_line = line.strip().split('\t', 1)
            code_set = line_map.get(rest_of_line, set())
            code_set.add(seq_code)
            line_map[rest_of_line] = code_set
        # write lines
        for line, code_set in line_map.items():
            out_file.write("{}\t{}\n".format(','.join(code_set), line))
        out_file.flush()


if __name__ == "__main__":
    os.chdir('D:\\Matan\\Dropbox\\PHd\\RiboSearch\\Mini project')
    rem_dups('match_log_nocm_14_4_18_no_dup', 'match_log')
