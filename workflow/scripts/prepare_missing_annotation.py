import sys
from Bio import SeqIO
import pandas as pd

##########################################################################

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

##########################################################################

tsv_df = pd.read_table(snakemake.input.tsv)

dict_translation = pd.read_table(
    snakemake.input.translation_table, index_col=1
).old_contig_name.to_dict()

all_wanted = tsv_df[(tsv_df.virsorter_cat.isna()) | (tsv_df.virsorter_cat == "No")]
all_wanted = all_wanted.contig_id.tolist()

with open(snakemake.output.fasta, "wt") as w_file:
    parser = SeqIO.parse(snakemake.input.fasta, "fasta")

    for contig in parser:
        # Because name in selected tsv and fasta file are different, the fasta file have {name_contig}-contigs- at the begining
        contig_id = dict_translation[contig.id]

        print(all_wanted)
        print(contig.id)
        print(contig_id)
        print(contig_id in all_wanted)

        # If prophage detected by checkv
        contig_id_prophage = "_".join(contig_id.split("_")[:-1])

        if contig_id in all_wanted or contig_id_prophage in all_wanted:
            SeqIO.write(contig, w_file, "fasta")
