from importlib_metadata import metadata
import pandas as pd
import sys
import os

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

# All files from the different blast databases
merge_blast_file = snakemake.input.merge_blast

# Metadata table if exists
metadata_file = snakemake.params.metadata

if os.path.isfile(metadata_file):
    metadata_table = pd.read_table(metadata_file)
    metadata_databases = metadata_table.database.unique().tolist()
else:
    metadata_table = ''
    metadata_databases = []

# Read blast file
big_blast = pd.read_table(merge_blast_file)

big_blast = big_blast.merge(full_meta_database.rename(columns={'contig_id': 'hit_genome'}), on=['hit_genome', 'database'], how='left')

new_blast = []
for index, g in big_blast.groupby('database'):
    columns2keep_g = ['contig', 'viral_taxonomy', 'host_taxonomy']
    g = g.sort_values(["evalue", "coverage", "pident"]).drop_duplicates(["contig"])
    g = g[columns2keep_g]
    g = g.rename(columns={i:f'{i}_{index}' for i in columns2keep_g[1:]})
    g = g.dropna(how='all', axis=1)

    if g.shape[1] == 1:
        g[f"present_in_{index}"] = "yes"
    
    new_blast.append(g.set_index('contig'))