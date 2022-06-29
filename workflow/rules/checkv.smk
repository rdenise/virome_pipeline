##########################################################################
##########################################################################
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule checkv_setup:
    output:
        directory(
            os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "checkv_db",
                "checkv-db-v1.2",
            )
        ),
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "checkv", "checkv_setup.log"),
    params:
        checkv_db=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "checkv_db",
        ),
    resources:
        cpus=1,
    conda:
        "../envs/checkv.yaml"
    threads: 1
    shell:
        """
        checkv download_database '{params.checkv_db}' &> '{log}'
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
            "vs2-step1",
            "final-viral-combined.fa",
        ),
        database=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "checkv_db",
            "checkv-db-v1.2",
        ),
    output:
        viruses=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "virsorter",
            "{sample}",
            "viruses.fna",
        ),
        proviruses=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "virsorter",
            "{sample}",
            "proviruses.fna",
        ),
        contamination=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "virsorter",
            "{sample}",
            "contamination.tsv",
        ),
    params:
        output_dir=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "virsorter",
            "{sample}",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "checkv",
            "virsorter",
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
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule checkv_run_deepvirfinder:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "deepvirfinder",
            "{sample}",
            "{sample}.extract.fasta",
        ),
        database=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "checkv_db",
            "checkv-db-v1.2",
        ),
    output:
        viruses=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "deepvirfinder",
            "{sample}",
            "viruses.fna",
        ),
        proviruses=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "deepvirfinder",
            "{sample}",
            "proviruses.fna",
        ),
        contamination=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "deepvirfinder",
            "{sample}",
            "contamination.tsv",
        ),
    params:
        output_dir=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "deepvirfinder",
            "{sample}",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "checkv",
            "deepvirfinder",
            "{sample}.checkv_run.log"),
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