##########################################################################
##########################################################################
##
##                                Library
##
##########################################################################
##########################################################################

import os, sys
import pandas as pd
import numpy as np
from snakemake.utils import validate

##########################################################################
##########################################################################
##
##                               Functions
##
##########################################################################
##########################################################################


def get_final_output(outdir, contigs_list):
    """
    Generate final output name
    """
    final_output = []

    final_output += (
        os.path.join(OUTPUT_FOLDER, "results", "taxonomic_annotation_contigs.tsv"),
    )

    final_output += expand(
        os.path.join(
            outdir,
            "processing_files",
            "vcontact2",
            "{sample}",
            "genome_by_genome_overview.csv",
        ),
        sample=contigs_list,
    )

    final_output += (
        os.path.join(
            outdir,
            "processing_files",
            "hmmer",
            f"significant_hit.full_{hmm_evalue_full:.0e}.dom_{hmm_evalue_dom}.domtblout.txt",
        ),
    )

    final_output += expand(
        os.path.join(
            outdir, "results", "dramv", "distill", "{sample}", "amg_summary.tsv"
        ),
        sample=contigs_list,
    )

    final_output += expand(
        os.path.join(
            outdir,
            "results",
            "dramv",
            "distill",
            "{sample}",
            "merge",
            "amg_summary.tsv",
        ),
        sample=contigs_list,
    )

    return final_output


##########################################################################


def create_folder(mypath):
    """
    Created the folder that I need to store my result if it doesn"t exist
    :param mypath: path where I want the folder (write at the end of the path)
    :type: string
    :return: Nothing
    """

    try:
        os.makedirs(mypath)
    except OSError:
        pass

    return


##########################################################################


def max_len_seq(file_fasta, ext_compress):
    max_len = 0
    if ext_compress == "tar.gz":
        import tarfile

        with tarfile.open(file_fasta, "r:gz") as tar:
            for tarinfo in tar:
                f = tar.extractfile(tarinfo.name)
                # To get the str instead of bytes str
                # Decode with proper coding, e.g. utf-8
                content = f.read().decode("utf-8", errors="ignore")
                # Split the long str into lines
                # Specify your line-sep: e.g. \n
                lines = content.split("\n")

                for line in lines:
                    if line.startswith(">"):
                        max_len = max(max_len, tmp_len)
                        tmp_len = 0
                    else:
                        tmp_len += len(line)
    elif ext_compress == "gz":
        import gzip

        with gzip.open(file_fasta, "rt") as r_file:
            for line in r_file:
                if line.startswith(">"):
                    max_len = max(max_len, tmp_len)
                    tmp_len = 0
                else:
                    tmp_len += len(line)

    elif ext_compress == "":
        with open(file_fasta, "rt") as r_file:
            for line in r_file:
                if line.startswith(">"):
                    max_len = max(max_len, tmp_len)
                    tmp_len = 0
                else:
                    tmp_len += len(line)

    print(max_len)
    print(tmp_len)

    return max_len


##########################################################################
##########################################################################
##
##                                Variables
##
##########################################################################
##########################################################################

# Validation of the config.yaml file
validate(config, schema="../schemas/config.schema.yaml")

# path to database sheet (TSV format, columns: database_name, path_db)
db_file = config["databases"]

# Validation of the database file
db_dtypes = {
    "database_name": "string",
    "database_filename": "string",
    "path_db": "string",
    "db_format": "string",
}

db_table = pd.read_table(db_file, dtype=db_dtypes)

validate(db_table, schema="../schemas/databases.schema.yaml")

##########################################################################
##########################################################################
##
##                        Core configuration
##
##########################################################################
##########################################################################

## Store some workflow metadata
config["__workflow_basedir__"] = workflow.basedir
config["__workflow_basedir_short__"] = os.path.basename(workflow.basedir)
config["__workflow_workdir__"] = os.getcwd()

if workflow.config_args:
    tmp_config_arg = '" '.join(workflow.config_args).replace("=", '="')
    config["__config_args__"] = f' -C {tmp_config_arg}"'
else:
    config["__config_args__"] = ""

with open(os.path.join(workflow.basedir, "../config/VERSION"), "rt") as version:
    url = "https://github.com/vdclab/sORTholog/releases/tag"
    config["__workflow_version__"] = version.readline()
    config["__workflow_version_link__"] = f"{url}/{config['__workflow_version__']}"


##########################################################################
##########################################################################
##
##                           Options
##
##########################################################################
##########################################################################

# Name your project
project_name = config["project_name"]

# Result folder
OUTPUT_FOLDER = os.path.join(config["output_folder"], project_name)
# Adding to config for report
config["__output_folder__"] = os.path.abspath(OUTPUT_FOLDER)

# Options for blastn
blast_evalue = config["default_blast_option"]["e_val"]
blast_coverage = config["default_blast_option"]["coverage"]
blast_pident = config["default_blast_option"]["pident"]

# Options for prokka
prokka_protein_db = config["default_prokka_option"]["protein_db"]
prokka_hmm_db = config["default_prokka_option"]["hmm_db"]
prokka_kingdom = config["default_prokka_option"]["kingdom"].capitalize()

# Option for DeepVirFinder
cutoff_deepvirfinder = config["default_deepvirfinder_option"]["cutoff_length"]

# Option for virsorter
cutoff_virsorter = config["default_virsorter_option"]["cutoff_length"]

# Option for DRAMv
cutoff_dramv = config["default_dramv_option"]["cutoff_length"]

# Options for hmmer
hmm_evalue_full = config["default_hmmer_option"]["e_val_full"]
hmm_evalue_dom = config["default_hmmer_option"]["e_val_dom"]

# Annotation table for taxonomy validation
if config["annotation_phages"]:

    # path to database sheet (TSV format, columns: "contig_id", "viral_taxonomy", "host_taxonomy", "database")
    phage_annotation_file = config["annotation_phages"]

    # Validation of the database file
    phage_annotation_dtypes = {
        "contig_id": "string",
        "viral_taxonomy": "string",
        "host_taxonomy": "string",
        "database": "string",
    }

    phage_annotation_table = pd.read_table(
        phage_annotation_file, dtype=phage_annotation_dtypes, sep="\t"
    )

    validate(phage_annotation_table, schema="../schemas/annotation_phages.schema.yaml")


DB_DICT = {"hmm": {}, "fasta": {}, "human": {}}

# Create a dictionary of the database file order by format
for index, row in db_table.iterrows():
    database_name = row.database_name
    DB_DICT[row.db_format.lower()][database_name] = {
        "path": row.path_db,
        "file": row.database_filename,
    }

# path to contigs sheet (TSV format, columns: contig_name, path_contig)
CONTIGS_FOLDER = config["contigs"]

if not config["contigs_ext"].startswith("."):
    CONTIGS_EXT = f".{config['contigs_ext']}"
else:
    CONTIGS_EXT = config["contigs_ext"]

# Get all the files int the contigs folder
(CONTIGS_FILES,) = glob_wildcards(
    os.path.join(CONTIGS_FOLDER, "{contigs_files}" + CONTIGS_EXT)
)
CONTIGS_DICT = {}

EXT_COMPRESS = ""

# Create a dictionary of the contigs files
for contig_file in CONTIGS_FILES:
    contig_name = contig_file.split(".")[0]
    contig_name_file = contig_file + CONTIGS_EXT

    # Test if the contigs are compressed to uncompress in case
    if "tar.gz" in CONTIGS_EXT:
        EXT_COMPRESS = "tar.gz"
        contig_name_file = contig_name_file.replace(".tar.gz", "")
    elif "gz" in CONTIGS_EXT:
        EXT_COMPRESS = ".gz"
        contig_name_file = contig_name_file.replace(".gz", "")

    MAX_LEN = max_len_seq(
        os.path.join(CONTIGS_FOLDER, contig_name_file), ext_compress=EXT_COMPRESS
    )

    max_cutoff = max(cutoff_virsorter, cutoff_deepvirfinder)
    # Remove from the analysis files that have contig too short:
    if MAX_LEN >= max(cutoff_virsorter, cutoff_deepvirfinder):
        CONTIGS_DICT[contig_name] = {
            "file": contig_name_file,
            "ext_compress": EXT_COMPRESS,
        }
    else:
        print(f"The file: {contig_name_file} does not have contigs above {max_cutoff}")
