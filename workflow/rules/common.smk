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


def get_final_output(outdir, contigs_list, datbase_list, blast_evalue):
    """
    Generate final output name
    """
    final_output = []

    final_output += expand(
        os.path.join(
            outdir,
            "processing_files",
            "blast",
            "contig_{contig}.evalue_{evalue:.0e}.{database}.blastn.outfmt6.txt"
        ),
        contig = contigs_list,
        database = datbase_list,
        evalue = [blast_evalue]
    )

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

for index, row in db_table.iterrows():
    database_name = row.database_name.split(".")[0]
    DB_DICT[row.db_format.lower()] = {
                       "path":row.path_db,
                       "file":row.database_name,
                       }

# path to contigs sheet (TSV format, columns: contig_name, path_contig)
contigs_file = config["contigs"]

# Validation of the contig file
contigs_dtypes = {
    "contig_name": "string",
    "path_contig": "string",
}

contigs_table = pd.read_table(contigs_file, dtype=contigs_dtypes)

validate(contigs_table, schema="../schemas/contigs.schema.yaml")

CONTIGS_DICT = {}

for index, row in db_table.iterrows():
    contig_name = row.contig_name.split(".")[0]
    CONTIGS_DICT[contig_name] = {
                "path":row.path_contig, 
                "file":row.contig_name,
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

blast_evalue = config['default_blast_option']['e_val']
