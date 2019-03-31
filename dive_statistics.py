#!/usr/bin/env python3
"""
Generates statistics and images for dive results.
F
"""
import os
import shutil
from subprocess import Popen
import logging


def generate_r2r(folder_path: str, sto_file_name: str, force: bool=False):
    pre_comp_sto = os.path.join(folder_path, "{}_cons.sto".format(sto_file_name[:-4]))
    meta_file_path = os.path.join(folder_path, "{}.meta".format(sto_file_name[:-4]))
    result_img_path = os.path.join(folder_path, "{}.svg".format(sto_file_name[:-4]))
    if os.path.exists(result_img_path) and not force:
        return result_img_path
    try:
        args = ["r2r", "--GSC-weighted-consensus", os.path.join(folder_path, sto_file_name), pre_comp_sto,
                "3", "0.97", "0.9", "0.75", "4", "0.97", "0.9", "0.75", "0.5", "0.1"]
        with Popen(args) as proc:
            proc.communicate()
        if not os.path.exists(pre_comp_sto) or os.stat(pre_comp_sto).st_size == 0:
            logging.error("Failed to compute consensus structure for {}".format(sto_file_name))
            return None
        with open(meta_file_path, 'w') as meta_file:
            meta_file.write(pre_comp_sto)
            meta_file.flush()
        args = ["r2r", "--disable-usage-warning", meta_file_path, result_img_path]
        with Popen(args) as proc:
            proc.communicate()
        if not os.path.exists(result_img_path) or os.stat(result_img_path).st_size == 0:
            if os.path.exists(result_img_path) and os.stat(result_img_path).st_size == 0:
                os.remove(result_img_path)
            logging.error("Failed to generate svg image for {}".format(result_img_path))
            return None
    finally:
        if os.path.exists(pre_comp_sto):
            os.remove(pre_comp_sto)
        if os.path.exists(meta_file_path):
            os.remove(meta_file_path)
    logging.info("Created r2r for {}: {}".format(sto_file_name, result_img_path))
    return result_img_path


def generate_cmv(folder_path: str, cm_file_name: str, sto_file_name: str,
                 model_type: str='detailed', force: bool = False):
    out_svg = os.path.join(folder_path, "{}_{}.svg".format(cm_file_name[:-3], model_type))
    if os.path.exists(out_svg) and not force:
        return out_svg
    args = ["conda", "run", "CMV", "-m", os.path.join(folder_path, cm_file_name),
            '-s', os.path.join(folder_path, sto_file_name), '-d', model_type, '-f', 'svg',
            '-p', folder_path]
    with Popen(args) as proc:
        proc.communicate()
    normal_out_svg = os.path.join(folder_path, "{}.svg".format(cm_file_name[6:-3], model_type))
    if os.path.exists(normal_out_svg):
        shutil.move(normal_out_svg, out_svg)
    else:
        return None
    return out_svg


def generate_visuals(folder_path: str):
    files = os.listdir(folder_path)
    cm_files = [file_name for file_name in files if file_name[-3:] == '.cm']
    for cm_name in cm_files:
        sto_name = "{}.sto".format(cm_name[:-3])
        generate_cmv(folder_path, cm_name, sto_name)
        generate_r2r(folder_path, sto_name)


def generate_filter(folder_path: str, cuttoff: int):
    result = []
    with open(os.path.join(folder_path, "FINAL_summary"), 'r') as sum_file:
        sum_file.readline()
        for line in sum_file:
            values = line.strip().split('\t')
            if len(values) < 7:
                continue
            if int(values[2]) > cuttoff:
                result.append(values[0].strip())
    with open(os.path.join(folder_path, "Filter_list_{}.txt".format(cuttoff)), 'w') as out_file:
        out_file.write("\n".join(result))
        out_file.flush()
    return result


if __name__ == '__main__':
    threshold = [5, 10, 25, 50]
    base_folder = "/opt/home/matan/Dropbox/PHd/RiboSearch/SandD_finals"
    for cut in threshold:
        left_list = generate_filter(base_folder, cut)
    generate_visuals(base_folder)
