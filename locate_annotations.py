import os
from typing import Tuple, List, Dict, Optional
import gffutils
from dive_statistics import get_distance, change_dups, OUT_RANGE, get_gff
from compare_rfam import break_id





def nearby_gene(gff_path: str, start: int, end: int, strand: str) -> List[Tuple[str, str, str]]:
    gff_name = os.path.splitext(gff_path)[0]
    gff_db = "{}.db".format(gff_name)
    keep_building = True
    while keep_building:
        keep_building = False
        if not os.path.exists(gff_db):
            change_dups(gff_path)
            try:
                db = gffutils.create_db(gff_path, dbfn=gff_db, from_string=False)
            except ValueError as ve:
                if 'empty file provided' in str(ve):
                    try:
                        os.remove(gff_db)
                    except:
                        pass
                    return []
                else:
                    print("failed to create {}".format(gff_db))
                    raise
            except:
                print("failed to create {}".format(gff_db))
                raise
        else:
            try:
                db = gffutils.FeatureDB(gff_db)
            except:
                os.remove(gff_db)
                keep_building = True
    results = []
    for feature in db.region(start=max(0, min(start, end) - OUT_RANGE), end=max(start, end) + OUT_RANGE,
                             strand=strand):
        new_distance, new_relation = get_distance(start, end, feature)
        gene_id = feature.attributes.get('ID')[0]
        desc = feature.attributes.get('product')
        if desc is None:
            desc = feature.attributes.get('mol_type')
        results.append((str(gene_id), str(feature.featuretype), str(new_distance), str(new_relation), str(desc)))
    if not results:
        for feature in db.region(start=max(0, min(start, end) - OUT_RANGE), end=max(start, end) + OUT_RANGE,
                                 strand='+'):
            new_distance, new_relation = get_distance(start, end, feature)
            gene_id = "REVERSE:" + str(feature.attributes.get('ID')[0])
            desc = feature.attributes.get('product')
            if desc is None:
                desc = feature.attributes.get('mol_type')
            results.append((str(gene_id), str(feature.featuretype), str(new_distance), str(new_relation), str(desc)))
    return results


def read_keys(folder_path: str, input_name: str) -> List[Tuple[str, int]]:
    results = []
    with open(os.path.join(folder_path, input_name), 'r') as input_file:
        input_file.readline()
        for line in input_file:
            if line.strip() != '':
                words = line.strip().split('\t')
                results.append((words[0], int(words[2])))
    return results


if __name__ == "__main__":
    BASE_DIR = "/opt/home/matan/Dropbox/PHd/RiboSearch/redocm"
    '''
    name = "transketolase_ext_new"
    in_name = "transketolase_ext_new.txt"
    name = "156_2_myosin"
    in_name = "156_2_myosin.txt"
    '''
    name = "M_6_mammalia"
    in_name = "M_6_mammalia.txt"
    with open(os.path.join(BASE_DIR, "{}_all_annotations.txt".format(name)), 'w') as out_file:
        out_file.write('Match ID\ttaxonomy\tannotation ID\type\tdistance\trelation\tproduct\n')
        for item_id, item_taxonomy in read_keys(BASE_DIR, in_name):
            item_dict = break_id(item_id)
            gff_path = get_gff(item_dict['code'])
            out_file.write("\n")
            for annotation in nearby_gene(gff_path, int(item_dict['start']),
                                          int(item_dict['end']), item_dict['strand']):
                out_file.write("{}\t{}\t{}\n".format(item_id, item_taxonomy, "\t".join(annotation)))
