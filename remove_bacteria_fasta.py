import result_dive
from Bio import SeqIO

def remove_seq_by_ancestor(fasta_in_path:str, fasta_out_path:str, ancestor_name:str='Bacteria'):
    with open(fasta_in_path, 'rU') as fasta_in, open(fasta_out_path, 'w') as fasta_out:
        for record in SeqIO.parse(fasta_in, "fasta"):
            org_id = result_dive.get_tax_id(record.id)
            if not result_dive.check_ancestor(ancestor_name, org_id):
                fasta_out.write(">{}\n{}\n".format(record.description, record.seq))
                #print(">{}\n{}\n".format(record.description, record.seq))


if __name__ == "__main__":
    #remove_seq_by_ancestor("test.fasta", "test_nobac.fasta")
    remove_seq_by_ancestor("/DB/fasta_db/nt/nt", "/DB/fasta_db/nt/nt_nobac")
