##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule dram_setup:
    output:
        os.path.join(
            "databases",
            "dram_db",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "dram",
            "dram_setup.log"
        ),
    resources:
        cpus=1,
    conda:
        "../envs/dram.yaml"
    threads: 1
    shell:
        """
        DRAM-setup.py prepare_databases --skip_uniref --output_dir '{output}' &> '{log}'
        """


##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule dramv_annotate:
    input:
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
        database=os.path.join(
            "databases",
            "dram_db",
        ),
    output:
        tsv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "annotations.tsv"
        ),  
    params:
        output_dir=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
        ),
        cutoff=cutoff_dramv,
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "dramv",
            "{sample}.dramv_annotate.log"
        ),
    resources:
        cpus=5,
    conda:
        "../envs/dram.yaml"
    threads: 5
    shell:
        """
        DRAM-v.py annotate -i '{input.fasta}' -v '{input.viral_affi}' \
        -o '{params.outpur_dir}' --skip_trnascan --threads {threads} \
        --min_contig_size '{parms.cutoff}' &> '{log}'
        """


##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule dramv_distill:
    input:
        tsv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "annotations.tsv"
        ), 
        database=os.path.join(
            "databases",
            "dram_db",
        ),
    output:
        amg_summary=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "distill",
            "{sample}",
            "amg_summary.tsv",
        ),
        viral_genome_summary=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "distill",
            "{sample}",
            "viral_genome_summary.tsv",
        ),
    params:
        output_dir=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "distill",
            "{sample}",
        ),
        cutoff=cutoff_dramv,
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "dramv",
            "{sample}.dramv_distill.log"
        ),
    resources:
        cpus=5,
    conda:
        "../envs/dram.yaml"
    threads: 5
    shell:
        """
        DRAM-v.py distill -i {input.tsv} -o '{params.output_dir}' &> '{log}'
        """


##########################################################################
##########################################################################

