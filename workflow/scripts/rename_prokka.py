import glob
import os
import re
from Bio import SeqIO

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

######################################################


def getTranslation(tbl, all_names, tsv_output):
    """get name translation for contig files from prokka

    Args:
        tbl (string): Path to the tbl file
        all_names (list of string): All the name in the fasta in order
        tsv_output (string): Path of the output tsv table with the prokka and contig name

    Returns:
        dict: Dictionnary that contains the prokka name in keys and the contig name in values
    """
    dictTrans = {}
    index = 0

    patternContigTag = re.compile(r"Prokka\|([A-Z]+_[0-9]+)")
    patternGeneTag = re.compile(r"locus_tag\t([A-Z]+_[0-9]+)\n")

    # Have a value to know if we expected locus_tag
    gene = False

    with open(tbl) as r_file:
        with open(tsv_output, "wt") as w_file:

            w_file.write("prokka_name\tname_in_contig\n")

            for line in r_file:
                if line.startswith(">Feature"):
                    locusTag = patternContigTag.search(line).group(1)

                    dictTrans[locusTag] = all_names[index]
                    w_file.write(f"{locusTag}\t{all_names[index]}\n")

                    index += 1
                    occurence = 1
                elif "gene" in line:
                    # Gene instance so set to true
                    gene = True
                elif gene and patternGeneTag.search(line):
                    oldname = patternGeneTag.search(line).group(1)
                    newname = f"{dictTrans[locusTag]}_{str(occurence).zfill(5)}"

                    occurence += 1

                    dictTrans[oldname] = newname
                    w_file.write(f"{oldname}\t{newname}\n")

                    # Find it so no need anymore
                    gene = False

    return dictTrans


######################################################

tbl_file = snakemake.input.tbl
folder_prokka = os.path.dirname(tbl_file)

all_files = glob.glob(tbl_file.replace("tbl", "*"))

fasta_contig = snakemake.input.fasta_contig

parser = SeqIO.parse(fasta_contig, "fasta")
allNames = [seq.id for seq in parser]

# Create a dictionary with the name to change
rep = getTranslation(
    tbl=tbl_file, all_names=allNames, tsv_output=snakemake.output.transTbl
)

# Add the escape to the regex sensitive character
rep = dict((re.escape(k), v) for k, v in rep.items())

# Create a string that contain all the possibilities
pattern = re.compile("|".join(rep.keys()))

for file2change in all_files:
    # Read the whole file
    with open(file2change) as r_file:
        oldstring = r_file.read()

    # Change all the names
    newstring = pattern.sub(lambda m: rep[re.escape(m.group(0))], oldstring)

    # Write in the same file
    with open(file2change.replace("prokka_name", "contig_renamed"), "wt") as w_file:
        w_file.write(newstring)

######################################################
