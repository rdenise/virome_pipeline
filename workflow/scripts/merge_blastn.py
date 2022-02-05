import sys
import os

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

blastn_files = snakemake.input.all_out

with open(snakemake.output.tsv, "w") as w_file:
    header = "qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle".replace(' ', '\t')
    w_file.write(f"{header}\n")

    for blast_file in blastn_files:
        with open(blast_file) as r_file:

            # Filter the blast hits using threshold 
            evalue = snakemake.wildcards.evalue
            coverage = snakemake.wildcards.coverage

            for line in r_file:
                line_split = line.split()

                evalue_blast = float(line_split[6])
                qcoverage = float(line_split[8])-float(line_split[7]) / float(line_split[4])
                #scoverage = float(line_split[10])-float(line_split[9]) / float(line_split[5])

                if evalue_blast <= evalue and coverage >= qcoverage :
                    w_file.write(line)


###########################################################
###########################################################
