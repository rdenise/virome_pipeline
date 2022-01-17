##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule virsorter_setup:
    output:
        os.path.join(
            "databases",
            "virsorter_db",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "virsorter",
            "virsorter_setup.log"
        ),
    resources:
        cpus=4,
    conda:
        "../envs/virsorter.yaml"
    threads: 4
    shell:
        """
        virsorter setup -d '{output}' -j '{threads}' &> '{log}'
        """


##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule virsorter_run:
    input:
        contig=lambda wildcards: os.path.join(
            CONTIGS_FOLDER,
            CONTIGS_DICT[wildcards.sample]["file"],
        ),
        database=os.path.join(
            "databases",
            "virsorter_db",
        ),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "final-viral-combined.fa",
        ),
        score=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "final-viral-score.tsv",
        ),
        boundary=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "final-viral-boundary.tsv",
        ),    
    params:
        output_dir=directory(os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
        )),           
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "virsorter",
            "{sample}.virsorter_run.log"
        ),
    resources:
        cpus=5,
    conda:
        "../envs/virsorter.yaml"
    threads: 5
    shell:
        """
        virsorter run -w '{params.output_dir}' -i '{input.contig}' -j {threads} all &> '{log}'
        """


##########################################################################
##########################################################################

