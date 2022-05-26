import pandas as pd
import sys
import os
from pathlib import Path

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

# All files from the different blast databases
merge_blast_file = snakemake.input.merge_blast

# Read blast file
big_blast = pd.read_table(merge_blast_file)

# Metadata table if exists
metadata_file = snakemake.params.metadata

if os.path.isfile(metadata_file):
    metadata_table = pd.read_table(metadata_file)

    # Add metadata
    big_blast["hit_genome"] = big_blast.hit_genome.replace(r"^(NC_[0-9]+)\.[0-9]$", r"\1", regex=True)
    
    big_blast = big_blast.merge(
        metadata_table.rename(columns={"contig_id": "hit_genome"}),
        on=["hit_genome", "database"],
        how="left",
    )

new_blast = []
for index, g in big_blast.groupby("database"):
    columns2keep_g = ["contig", "viral_taxonomy", "host_taxonomy"]
    g = g.sort_values(["evalue", "coverage", "pident"]).drop_duplicates(["contig"])
    g = g[columns2keep_g]
    g = g.rename(columns={i: f"{i}_{index}" for i in columns2keep_g[1:]})
    g = g.dropna(how="all", axis=1)

    if g.shape[1] == 1:
        g[f"present_in_{index}"] = "yes"

    new_blast.append(g.set_index("contig"))

all_contigs = pd.concat(new_blast, axis=1).reset_index().sort_values("contig")

# Load viral detection informations
path_tsv = Path(snakemake.params.viral_tsv)

tsv_detection = []

for tsv in path_tsv.glob("*.selected.tsv"):
    tsv_df = pd.read_table(tsv).rename(columns={"contig_id": "contig"})
    tsv_df["contig"] = tsv_df.contig.apply(
        lambda x: f"{str(tsv.stem).replace('.selected', '')}-{x}".replace("_", "-")
    )
    tsv_detection.append(tsv_df)

tsv_detection = pd.concat(tsv_detection)

# The most complete names or in the right file
all_contigs = (
    all_contigs.merge(tsv_detection, on="contig", how="right")
    .fillna("")
    .sort_values("contig", ignore_index=True)
)

# Load additional information from virsorter
path_virsorter = Path(snakemake.params.virsorter_tsv)

tsv_virsorter = []
wanted_columns = [
    "contig",
    "vs2_completeness",
    "num_hallmark_genes",
    "perc_viral_genes",
    "perc_non_viral_genes",
]

for tsv in path_virsorter.glob("*/vs2-step2/final-viral-score.tsv"):
    tsv_df = pd.read_table(tsv).rename(
        columns={
            "seqname": "contig",
            "viral": "perc_viral_genes",
            "cellular": "perc_non_viral_genes",
            "hallmark": "num_hallmark_genes",
        }
    )
    tsv_df["vs2_completeness"] = tsv_df.contig.apply(lambda x: x.split("|")[-1])
    tsv_df["contig"] = tsv_df.contig.apply(
        lambda x: f"{str(tsv.parts[-3])}-{x.split('|')[0]}".replace("_", "-")
    )
    tsv_df = tsv_df[wanted_columns]
    tsv_virsorter.append(tsv_df)

tsv_virsorter = pd.concat(tsv_virsorter)

# Add this information to all informations
all_contigs = all_contigs.merge(tsv_virsorter, on="contig", how="left")

all_contigs.to_csv(snakemake.output.tsv, sep="\t", index=False)
