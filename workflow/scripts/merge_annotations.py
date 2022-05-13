import sys
from Bio import SeqIO
import pandas as pd

##########################################################################

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

##########################################################################

tsv_missing = snakemake.input.tsv_missing

faa_missing = snakemake.input.faa_missing
        
tsv_virsorter = snakemake.input.tsv_virsorter

faa_virsorter = snakemake.input.faa_virsorter

transl_table = snakemake.input.translation_table
transl_dict = pd.read_table(transl_table, index_col=0).new_contig_name.to_dict()

##########################################################################

missing_df = pd.read_table(tsv_missing, index_col=0)
virsorter_df = pd.read_table(tsv_virsorter, index_col=0)

concat_df = pd.concat([virsorter_df, missing_df])
concat_df.to_csv(snakemake.output.tsv, sep='\t')

##########################################################################

with open(snakemake.output.fasta, "wt") as w_file:
    for faa_file in [faa_virsorter, faa_missing]:
        parser = SeqIO.parse(faa_file, "fasta")
        for protein in parser:
            protein_id_split = protein.id.split('__')

            if len(protein_id_split) > 1:
                if protein_id_split[0] in transl_dict:
                    print(f"old name: {protein.id}")
                    protein.id = transl_dict[protein_id_split[0]] + '_' + protein_id_split[-1].split('_')[-1]
                    print(f"new name {protein.id}")
                    print("---------------------")
                else :
                    continue

            protein.name = protein.description = ''
            SeqIO.write(protein, w_file, "fasta")

