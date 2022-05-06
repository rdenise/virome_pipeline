import sys
from Bio import SeqIO
import pandas as pd

##########################################################################

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

##########################################################################

tsv_df = pd.read_table(snakemake.input.tsv)

all_wanted = tsv_df[(tsv_df.virsorter_cat.isna()) | ~(tsv_df.virsorter_cat.str.contains('keep'))]
all_wanted = all_wanted.contig_id.tolist()

with open(snakemake.output.fasta, 'wt') as w_file:
    parser = SeqIO.parse(snakemake.input.fasta, 'fasta')

    for contig in parser:
        # Because name in selected tsv and fasta file are different, the fasta file have {name_contig}-contigs- at the begining
        contig_id = contig.id.split('contigs-')[-1]

        if contig.id in all_wanted:
            SeqIO.write(contig, w_file, "fasta")

            