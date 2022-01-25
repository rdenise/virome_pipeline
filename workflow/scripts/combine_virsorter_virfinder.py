from Bio import SeqIO
import sys
import os

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

# List that will contains all the contigs to filter
all_contig_ids = []

# Get all the names from the virsorter keep2 list
ids_virsorter_keep2 = snakemake.input.ids_virsorter_keep2

with open(ids_virsorter_keep2) as r_file:
    r_file.readline()

    for line in r_file:
        all_contig_ids.append(line.rsplit())

# Get all the names from the virsorter keep1 list and remove redondant name
ids_virsorter_keep1 = snakemake.input.ids_virsorter_keep1

with open(ids_virsorter_keep1) as r_file:
    r_file.readline()

    for line in r_file:
        rsplit_line = line.rsplit()
        
        if rsplit_line not in all_contig_ids:
            all_contig_ids.append(rsplit_line)

# Get all the names from the deepvirfinder list and remove redondant name
ids_virfinder = snakemake.input.ids_virfinder

with open(ids_virfinder) as r_file:
    r_file.readline()

    for line in r_file:
        rsplit_line = line.rsplit()
        
        if rsplit_line not in all_contig_ids:
            all_contig_ids.append(rsplit_line)

# Parse the fasta of the contig and create the new one
fasta_contigs = snakemake.input.contigs

with open(snakemake.output.fasta, "w") as w_file:

    parser = SeqIO.parse(fasta_contigs, "fasta")

    for contig in parser:
        if contig.id in all_contig_ids:
            contig.id = f"{snakemake.wildcards.sample}-{contig.id}".replace("_", "-")
            contig.name = ""
            config.description = ""

            SeqIO.write(contig, w_file, "fasta")
            


###########################################################
###########################################################
