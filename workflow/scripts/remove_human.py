import sys
import os
from Bio import SeqIO

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

blastn_file = snakemake.input.tsv

with open(snakemake.params.tsv, "w") as w_file:
    header = "qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle".replace(' ', '\t')
    w_file.write(f"{header}\n")

    contig2remove = []

    with open(blastn_file) as r_file:

        # Filter the blast hits using threshold 
        split_name = snakemake.params.tsv.split('.human')[0].split(snakemake.wildcards.sample)[-1]

        evalue = float(split_name.split('.cov')[0].split('.eval_')[-1])
        coverage = float(split_name.split('.cov_')[-1])

        for line in r_file:
            line_split = line.split()

            evalue_blast = float(line_split[6])
            qcoverage = (float(line_split[8])-float(line_split[7]) +_1) / float(line_split[4])
            #scoverage = (float(line_split[10])-float(line_split[9]) + 1) / float(line_split[5])

            if evalue_blast <= evalue and coverage >= qcoverage :
                w_file.write(line)
                contig2remove.append(line_split[0])

parser = SeqIO.parse(snakemake.input.fasta, 'fasta')

with open(snakemake.output.fasta, 'wt') as w_file:
    for contig in parser:
        if contig.id not in contig2remove:
            # To make sure that it is only a contig name as in the input fasta
            contig.name = contig.description = ''

            SeqIO.write(contig, w_file, 'fasta')

###########################################################
###########################################################
