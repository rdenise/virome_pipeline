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

    final_output += os.path.join(
            outdir,
            "processing_files",
            "blast",
            f"merge.eval_{blast_evalue:.0e}.cov_{blast_coverage}.annotation.blasn.tsv"
        ),

    final_output += expand(
        os.path.join(
            outdir,
            "processing_files",
            "vcontact2",
            "{sample}",
            "genome_by_genome_overview.csv",
        ),
        sample = contigs_list,
    )

    final_output += os.path.join(
            outdir,
            "processing_files",
            "hmmer",
            f"significant_hit.full_{hmm_evalue_full}.dom_{hmm_evalue_dom}.domtblout.txt"
        ),  

    return final_output


##########################################################################


def create_folder(mypath):
    """
    Created the folder that I need to store my result if it doesn't exist
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
    "path_db": "string",
    "db_format": "string",
}

db_table = pd.read_table(db_file, dtype=db_dtypes)

validate(db_table, schema="../schemas/databases.schema.yaml")

DB_DICT = {'hmm':{}, 'fasta':{}}

# Create a dictionary of the database file order by format
for index, row in db_table.iterrows():
    database_name = row.database_name.split(".")[0]
    DB_DICT[row.db_format.lower()][database_name] = {
                       "path":row.path_db,
                       "file":row.database_name,
                       }

# path to contigs sheet (TSV format, columns: contig_name, path_contig)
CONTIGS_FOLDER = config["contigs"]

if not config["contigs_ext"].startswith("."): 
    CONTIGS_EXT = f'.{config["contigs_ext"]}'
else: 
    CONTIGS_EXT = config["contigs_ext"]

# Get all the files int the contigs folder
CONTIGS_FILES, = glob_wildcards(os.path.join(
                                CONTIGS_FOLDER,
                                "{contigs_files}" +\
                                CONTIGS_EXT))

CONTIGS_DICT = {}

EXT_COMPRESS = ""

# Create a dictionary of the contigs files
for contig_file in CONTIGS_FILES:
    contig_name = contig_file.split(".")[0]
    contig_name_file = contig_file + CONTIGS_EXT

    # Test if the contigs are compressed to uncompress in case
    if "tar.gz" in CONTIGS_EXT:
        EXT_COMPRESS = "tar.gz"
        contig_name_file = contig_name_file.replace(
                    ".tar.gz", "")
    elif "gz" in CONTIGS_EXT:
        EXT_COMPRESS = ".gz"
        contig_name_file = contig_name_file.replace(
                    ".gz", "")

    CONTIGS_DICT[contig_name] = {
                "file":contig_name_file,
                "ext_compress":EXT_COMPRESS
        }

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
blast_evalue = config['default_blast_option']['e_val']
blast_coverage = config['default_blast_option']['coverage']

# Options for prokka
prokka_protein_db = config['default_prokka_option']['protein_db']
prokka_hmm_db = config['default_prokka_option']['hmm_db']
prokka_kingdom = config['default_prokka_option']['kingdom'].capitalize()

# Option for DeepVirFinder
cutoff_deepvirfinder = config['default_deepvirfinder_option']['cutoff_length']

# Option for virsorter
cutoff_virsorter = config['default_virsorter_option']['cutoff_length']

# Option for DRAMv
cutoff_dramv = config['default_dramv_option']['cutoff_length']

# Options for hmmer
hmm_evalue_full = config['default_hmmer_option']['e_val_full']
hmm_evalue_dom = config['default_hmmer_option']['e_val_dom']


