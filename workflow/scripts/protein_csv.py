from Bio import SeqIO
import sys

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

protein_file = snakemake.input.proteins_fasta

parser = SeqIO.parse(protein_file, 'fasta')

contig_parser = SeqIO.parse(snakemake.params.fasta_contig, "fasta")

dict_contig_len = {contig.id:len(contig.seq) for contig in contig_parser}

with open(snakemake.output.csv, 'wt') as w_file:
    with open(snakemake.output.fasta, 'wt') as fasta_file:
        with open(snakemake.output.fasta_low, 'wt') as fasta_file_low:
            header = "contig_id,protein_id,keywords"
            w_file.write(f"{header}\n")

            for protein in parser:
                protein_id = protein.id
                contig_id = protein_id.split("_")[0]
                keyword = 'None'

                if dict_contig_len[contig_id] >= 10000:
                    w_file.write(f"{contig_id},{protein_id},{keyword}\n")

                    protein.name = ""
                    protein.description = ""

                    SeqIO.write(protein, fasta_file, "fasta")
                else :
                    protein.name = ""
                    protein.description = ""

                    SeqIO.write(protein, fasta_file_low, "fasta")                

###########################################################
###########################################################
