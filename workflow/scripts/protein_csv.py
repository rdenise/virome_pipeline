# To convert in python

echo "contig_id,protein_id,keywords" > proteins.csv
grep '>' translations_vir_contigs1000.faa | sed 's/>//g' | sed 's/ # .*//g' > headers.tmp
cat headers.tmp | sed 's/$/,None/' > proteins.tmp
cat headers.tmp | sed 's/_[^_]*$/,/' > contigs.tmp
paste contigs.tmp proteins.tmp | sed 's/\t//g' >> proteins.csv
rm *tmp;
sed -i 's/ # .*//g' translations_vir_contigs1000.faa;

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
            protein_id = 
            contig_id = 
            keyword =

            w_file.write("{contig_id},{protein_id},None\n")

            protein.id =
            protein.name = ""
            protein.description = ""

            SeqIO.write(protein, fasta_file, "fasta")

###########################################################
###########################################################

pr