##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule virsorter_setup:
    output:
        directory(os.path.join(
            OUTPUT_FOLDER,
            "database",
            "virsorter_db"
        )),
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
            CONTIGS_DICT[wildcards.contig]["file"],
        ),
        database=directory(os.path.join(
            OUTPUT_FOLDER,
            "database",
            "virsorter_db"
        )),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "contig_{contig}",
            "final-viral-combined.fa",
        ),
        score=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "contig_{contig}",
            "final-viral-score.tsv",
        ),
        boundary=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "contig_{contig}",
            "final-viral-boundary.tsv",
        ),    
    params:
        output_dir=directory(os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "contig_{contig}",
        )),           
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "virsorter",
            "contig_{contig}.virsorter_run.log"
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
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule virsorter_postprocess:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "contig_{contig}",
            "final-viral-combined.fa",
        ),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "contig_{contig}",
            "virsorter_positive.ids",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "virsorter",
            "contig_{contig}.virsorter_postprocess.log"
        ),
    resources:
        cpus=1,
    threads: 1
    shell:
        """
        grep '>' '{input}' | cut -d '|' -f 1 | tr -d '>' > '{output}'
        """


##########################################################################
##########################################################################
