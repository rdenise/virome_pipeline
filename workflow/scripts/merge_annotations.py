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

##########################################################################

missing_df = pd.read_table(tsv_missing, index_col=0)
virsorter_df = pd.read_table(tsv_virsorter, index_col=0)

concat_df = pd.concat([virsorter_df, missing_df])
concat_df.to_csv(snakemake.output.tsv)

##########################################################################

with open(snakemake.output.fasta, "wt") as w_file:
    for faa_file in [faa_virsorter, faa_missing]:
        parser = SeqIO.parse(faa_file, "fasta")
        for protein in parser:
            protein.name = ''
            SeqIO.write(protein, w_file, "fasta")
