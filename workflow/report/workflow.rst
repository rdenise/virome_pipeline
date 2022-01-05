============================
Gene sORTholog report
============================

This virome_pipeline_ is a Snakemake workflow that ???. This analysis is based on commit version {{ snakemake.config["__workflow_version__"] }}_.

The analysis can be rerun with the following command:

Local:

.. code-block:: shell

   cd {{ snakemake.config["__workflow_workdir__"] }}
{% if snakemake.config["__workflow_basedir__"] == snakemake.config["__workflow_workdir__"] %}
   snakemake -j 1 --use-conda {{ snakemake.config["__config_args__"] }}
{% else %}
   snakemake -j 1 --use-conda {{ snakemake.config["__config_args__"] }} -s {{ snakemake.config["__workflow_basedir_short__"] }}/Snakefile
{% endif %}

.. note::

   Since the workflow is still work in progress, make sure to first 
   run the commands with the `--dry-run` (`-n`) flag to make sure you 
   don't inadvertedly have to regenerate large parts of the results.
   Many workflow dependencies are complex and in particular when 
   running smaller parts of the workflow, unexpected things may 
   happen.  


Workflow summary
----------------

The workflow runs the following steps:

1. ???


Data organization
-----------------

.. code-block:: text

   {{ snakemake.config["project_name"] }}/                          <- top-level project folder
   │
   │
   ├── report.html                           <- This report file      
   │
   ├── logs                                  <- Collection of log outputs, e.g. from cluster managers
   │
   ├── databases                             <- Generated analysis database related files
   │   ├── all_taxid                         
   │   │   ├─ protein_table.tsv              <- Table with the informations about the proteins of the downloaded taxid
   │   │   ├─ summary_assembly_taxid.tsv     <- Table with the informations about the downloaded genome from NCBI
   │   │   ├─ taxid_all_together.fasta       <- Fasta file of the downloaded taxid
   │   │   └─ taxid_checked.txt              <- List of the downloaded taxid
   │   │
   │   ├── merge_fasta                       
   │   │   └─ taxid_checked.txt              <- Fasta with the concatenation of the genome of interest and seeds
   │   │
   │   └── seeds                             
   │       ├─ seeds.fasta                    <- Fasta file of the seeds
   │       └─ new_seeds.tsv                  <- Table with the informations about the seeds
   │
   └── results                               <- Final results for sharing with collaborators, typically derived from analysis sets
       ├── patab_melt.tsv                    <- Table with the information of sORTholog one information by line
       ├── patab_table.tsv                   <- Table with the information of presence absence with genome in index and seeds in columns and proteins Id in the cell
       └── plots                             <- Plots and table on which the plot are created



Analysis overview
-----------------

The analyses can basically be divided in two parts: `Raw data
analysis`_ and `Analysis sets`_.

Raw data analysis
*****************

- In the first step it perform ???

Analysis sets
*************

???

General results
---------------

????

Workflow graph
--------------


.. _virome_pipeline: https://github.com/vdclab/virome_pipeline
.. _{{ snakemake.config["__workflow_version__"] }}: {{ snakemake.config["__workflow_version_link__"] }}
