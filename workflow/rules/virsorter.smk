##########################################################################
##########################################################################
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule virsorter_setup:
    output:
        directory(
            os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "virsorter_db",
            )
        ),
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "virsorter", "virsorter_setup.log"),
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


rule virsorter_run_step1:
    input:
        contig=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "contigs",
            "human_filtered",
            "{sample}.filtered.sorted.fasta",
        ),
        database=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "virsorter_db",
        ),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "vs2-step1",
            "final-viral-combined.fa",
        ),
        score=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "vs2-step1",
            "final-viral-score.tsv",
        ),
        boundary=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "vs2-step1",
            "final-viral-boundary.tsv",
        ),
    params:
        output_dir=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "vs2-step1",
        ),
        cutoff=cutoff_virsorter,
    log:
        os.path.join(
            OUTPUT_FOLDER, "logs", "virsorter", "{sample}.virsorter_run_step1.log"
        ),
    resources:
        cpus=5,
    conda:
        "../envs/virsorter.yaml"
    threads: 5
    shell:
        """
        virsorter run --keep-original-seq -i '{input.contig}' -w '{params.output_dir}' \
        --include-groups "dsDNAphage,ssDNA" --min-length {params.cutoff} \
        --min-score 0.9 -j {threads} all &> '{log}'
        """


##########################################################################
##########################################################################
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule virsorter_run_step2:
    input:
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
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "vs2-step2",
            "for-dramv",
            "final-viral-combined-for-dramv.fa",
        ),
        viral_affi=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "vs2-step2",
            "for-dramv",
            "viral-affi-contigs-for-dramv.tab",
        ),
    params:
        output_dir=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "vs2-step2",
        ),
        cutoff=cutoff_virsorter,
        input_vs2=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "virsorter",
            "{sample}",
            "combined.fna",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER, "logs", "virsorter", "{sample}.virsorter_run_step2.log"
        ),
    resources:
        cpus=5,
    conda:
        "../envs/virsorter.yaml"
    threads: 5
    shell:
        """
        cat '{input.viruses}' '{input.proviruses}' > '{params.input_vs2}'

        virsorter run --seqname-suffix-off --viral-gene-enrich-off \
        --provirus-off --prep-for-dramv -i '{params.input_vs2}' \
        -w '{params.output_dir}' --include-groups "dsDNAphage,ssDNA" \
        --min-length {params.cutoff} --min-score 0.5 \
        -j {threads} all &> '{log}'
        """


##########################################################################
##########################################################################
