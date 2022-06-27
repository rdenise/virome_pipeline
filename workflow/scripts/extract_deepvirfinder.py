from Bio import SeqIO
import sys

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

# Get all the names from the deepvirfinder list and remove redondant name
ids_virfinder = snakemake.input.ids_virfinder


# Parse the fasta of the contig and create the new one
fasta_contigs = snakemake.input.contig


with open(snakemake.output.fasta, "w") as w_file:
    with open(ids_virfinder) as r_file:
        r_file.readline()

        index_fasta = SeqIO.index(fasta_contigs, "fasta")

        for line in r_file:
            SeqIO.write(index_fasta[line.rstrip()], w_file, "fasta")

###########################################################
###########################################################
