##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule checkv_setup:
    output:
        os.path.join(
            "databases",
            "checkv_db",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "checkv",
            "checkv_setup.log"
        ),
    resources:
        cpus=1,
    conda:
        "../envs/checkv.yaml"
    threads: 1
    shell:
        """
        checkv download_database '{output}' &> '{log}'
        """


##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule checkv_run:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "final-viral-combined.fa",
        ),
        database=os.path.join(
            "databases",
            "checkv_db",
        ),
    output:
        viruses=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "{sample}",
            "viruses.fna",
        ),
        proviruses=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "{sample}",
            "proviruses.fna",
        ),    
        contamination=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "{sample}",
            "contamination.tsv",
        ),
        params:
        output_dir=directory(os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "{sample}",
        )),           
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "checkv",
            "{sample}.checkv_run.log"
        ),
    resources:
        cpus=5,
    conda:
        "../envs/checkv.yaml"
    threads: 5
    shell:
        """
        checkv end_to_end '{input.fasta}' '{params.output_dir}' \
        -t {threads} -d '{input.database}' &> '{log}'
        """


##########################################################################
##########################################################################

