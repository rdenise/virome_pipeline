import pandas as pd
import sys
import os

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################
# From the method for virsorter 2:
# Keep1: viral_gene >0
# Keep2: viral_gene =0 AND (host_gene =0 OR score >=0.95 OR hallmark >2)
# Manual check: (NOT in Keep1 OR Keep2) AND viral_gene =0 AND host_gene =1 AND length >=10kb
# Discard: the rest

# To look at the viral_gene, host_gene, score, and hallmark of sequences you can merge "vs2-pass1/final-viral-score.tsv" and "checkv/contamination.tsv", and filter in a spreadsheet.

# Get the informations from virsorter
deepvirfinder_file = snakemake.input.txt
deepvirfinder_file_df = pd.read_table(deepvirfinder_file)

# Good hit virfinder (form Camarillo-Gurerro_2021_Cell): score > 0.9 and pvalue < 0.01
deepvirfinder_file_df[
    (deepvirfinder_file_df.score > 0.9) & (deepvirfinder_file_df.pvalue < 0.01)
].name.to_csv(snakemake.output.txt, index=False)

###########################################################
###########################################################
