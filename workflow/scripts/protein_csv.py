from Bio import SeqIO
import sys
import os

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

protein_file = snakemake.input.proteins_fasta

parser = SeqIO.parse(protein_file, 'fasta')

with open(snakemake.output.csv, 'wt') as w_file:
    with open(snakemake.output.fasta, 'wt') as fasta_file:
        header = "contig_id,protein_id,keywords"
        w_file.write(f"{headet}\n")

        for protein in parser:
            # TO DO Because don't know the format and could not infer it form sed commands
            # Will be needed to modify the fasta file too because need to match the table probably
            protein_id = protein.id
            contig_id = snakemake.wildcards.contig
            keyword = 'None'

            w_file.write("{contig_id},{protein_id},{keyword}\n")

            protein.name = ""
            protein.description = ""

            SeqIO.write(protein, fasta_file, "fasta")

###########################################################
###########################################################

pr