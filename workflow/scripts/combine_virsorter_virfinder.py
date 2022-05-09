from Bio import SeqIO
import pandas as pd
import sys
import os

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

# List that will contains all the contigs to filter
all_contig_ids = []

# Dataframe that contains all the informations about
output_df = pd.DataFrame(columns=['contig_id', 'virsorter_cat', 'deepvirfinder'])

# Get all the names from the virsorter keep2 list
ids_virsorter_keep2 = snakemake.input.ids_virsorter_keep2_checked

with open(ids_virsorter_keep2) as r_file:
    r_file.readline()

    for line in r_file:
        rstrip_line = line.rstrip()
        rstrip_line = rstrip_line.split('||')[0]

        all_contig_ids.append(rstrip_line)

        output_df.at[rstrip_line, "contig_id"] = rstrip_line
        output_df.at[rstrip_line, "virsorter_cat"] = "keep2_checked"

# Get all the names from the virsorter keep1 list and remove redondant name
ids_virsorter_keep1 = snakemake.input.ids_virsorter_keep1

with open(ids_virsorter_keep1) as r_file:
    r_file.readline()

    for line in r_file:
        rstrip_line = line.rstrip()
        rstrip_line = rstrip_line.split('||')[0]
        
        if rstrip_line not in all_contig_ids:
            all_contig_ids.append(rstrip_line)

            output_df.at[rstrip_line, "contig_id"] = rstrip_line
            output_df.at[rstrip_line, "virsorter_cat"] = "keep1"

# Get all the names from the deepvirfinder list and remove redondant name
ids_virfinder = snakemake.input.ids_virfinder

with open(ids_virfinder) as r_file:
    r_file.readline()

    for line in r_file:
        rstrip_line = line.rstrip()

        output_df.at[rstrip_line, "contig_id"] = rstrip_line
        output_df.at[rstrip_line, "deepvirfinder"] = "Yes"

        if rstrip_line not in all_contig_ids:
            all_contig_ids.append(rstrip_line)

# Fill the informations missing now the list of contigs we keep is set 
dict_map_virsorter = {}

files_with_info = {snakemake.input.ids_virsorter_keep2_suspicious:"keep2_suspicious",
                   snakemake.input.ids_virsorter_manual_check:"to_manual_check",
                   snakemake.input.ids_virsorter_discarded:"discarded",
                   }

for file_ids in files_with_info:
    with open(file_ids) as r_file:
        r_file.readline()

        for line in r_file:
            rstrip_line = line.rstrip()
            rstrip_line = rstrip_line.split('||')[0]
            
            if rstrip_line not in all_contig_ids:
                dict_map_virsorter[rstrip_line] = files_with_info[file_ids]

# Fill the dataframe
list_contig2add_virsorter_cat = list(dict_map_virsorter.keys())
output_df.loc[output_df.contig_id.isin(list_contig2add_virsorter_cat),
             'virsorter_cat'] = output_df.loc[output_df.contig_id.isin(list_contig2add_virsorter_cat), 'contig_id'].map(dict_map_virsorter)

output_df.fillna('No', inplace=True)

# Parse the fasta of the contig and create the new one
fasta_contigs = snakemake.input.contigs

with open(snakemake.output.fasta, "w") as w_file:
    with open(snakemake.output.translation_table, "w") as tsv_file:
        tsv_file.write("old_contig_name\tnew_contig_name\n")

        parser = SeqIO.parse(fasta_contigs, "fasta")

        for contig in parser:
            if contig.id in all_contig_ids:
                contig_id = f"{snakemake.wildcards.sample}-{contig.id}".replace("_", "-")

                tsv_file.write(f"{contig.id}\t{contig_id}\n")

                contig.id = contig_id
                contig.name = ""
                contig.description = ""

                SeqIO.write(contig, w_file, "fasta")
            
output_df.to_csv(snakemake.output.tsv, sep="\t", index=False)

###########################################################
###########################################################
