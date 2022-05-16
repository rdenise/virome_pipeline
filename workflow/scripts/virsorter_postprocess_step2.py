import pandas as pd
import sys
import os

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################
# From the method for virsorter 2:
# Keep2: viral_gene =0 AND (host_gene =0 OR score >=0.95 OR hallmark >2)
# subset the DRAMv table using contigs in the "Keep2" category, and screen for the "suspicious" genes in the subset DRAMv table

# To look at the annotation.tsv from dramv

# Get the informations from virsorter
virsorter_keep2 = snakemake.input.ids_keep2
keep2_names = []

with open(virsorter_keep2) as r_file:
    header = r_file.readline()

    for name in r_file:
        keep2_names.append(name.rstrip())

# Get the annotation form DRAMv
dramv_annot = snakemake.input.annotations
dramv_df = pd.read_table(dramv_annot, dtype="string")
dramv_df.fillna("", inplace=True)

# Change name dramv to virsorter
dramv_df["contig_id"] = dramv_df.scaffold.apply(
    lambda x: x.split("-cat")[0].replace("__", "||")
)

# Only keep the names of interest
dramv_df = dramv_df[dramv_df.contig_id.isin(keep2_names)].reset_index(drop=True)

# Suspicious genes
suspicious_gene = snakemake.input.suspicous_gene

suspicious_names = []

with open(suspicious_gene) as r_file:
    for name in r_file:
        suspicious_names.append(name.rstrip())

# Check for each columns if suspicious genes name exists
suspicious_index = dramv_df.fasta.str.contains(suspicious_names[0])

# For each suspicious gene for each columns in the dataframe search for suspicious genes
for name in suspicious_names:
    for column in dramv_df.columns:
        suspicious_index += dramv_df[column].str.contains(name)

# Get the name of all suspicious scaffolds
suspicious_scaffolds = dramv_df.loc[suspicious_index, "contig_id"].tolist()

# Reduction of the dataframe to be able to have only one scaffold name per line
dramv_df = dramv_df.drop_duplicates("contig_id").reset_index(drop=True)

# Write the suspicious in a file
dramv_df.loc[dramv_df.contig_id.isin(suspicious_scaffolds), "contig_id"].to_csv(
    snakemake.output.suspicious, sep="\t", index=False
)

# Write the checked in another file
dramv_df.loc[~(dramv_df.contig_id.isin(suspicious_scaffolds)), "contig_id"].to_csv(
    snakemake.output.checked, sep="\t", index=False
)

###########################################################
###########################################################
