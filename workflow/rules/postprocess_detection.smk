##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule deepvirfinder_postprocess:
    input:
        txt=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "deepvirfinder",
            "{sample}",
            "{sample}_gt{cutoff}bp_dvfpred.txt",
        ),
    output:
        txt=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "deepvirfinder",
            "{sample}",
            "{sample}.deepvirfinder_positive.gt{cutoff}bp.ids",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "deepvirfinder",
            "{sample}.{cutoff}.deepvirfinder_postprocess.log"
        ),
    resources:
        cpus=1,
    threads: 1
    script: 
        "../scripts/deepvirvinder_posprocess.py"


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
            "{sample}",
            "final-viral-combined.fa",
        ),
    output:
        ids=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter_positive.ids",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "virsorter",
            "{sample}.virsorter_postprocess.log",
        ),
    resources:
        cpus=1,
    conda:
        "../envs/biopython.yaml"
    threads: 1
    script: 
        "../scripts/virsorter_posprocess.py"


##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule combine_virsorter_virfinder:
    input:
        ids_virsoter=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter_positive.ids",
        ),
        ids_virfinder=lambda wildcards: os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "deepvirfinder",
            wildcards.sample,
            f"{wildcards.sample}.deepvirfinder_positive.gt{cutoff_deepvirfinder}bp.ids",
        ),        
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "viral_contigs",
            "{sample}.selected.fasta",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "postprocess_detection",
            "{sample}.combine_virsorter_virfinder.log"
        ),
    resources:
        cpus=1,
    conda:
        "../envs/biopython.yaml"
    threads: 1
    script: 
        "../scripts/combine_virsorter_virfinder.py"


##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule merge_blastn:
    input:
        all_out=expand(
            os.path.join(
                OUTPUT_FOLDER,
                "processing_files",
                "blast",
                "{sample}.evalue_{evalue:.0e}.{database}.blastn.outfmt6.txt"
            ),
            sample = CONTIGS_DICT.keys(),
            database = DB_DICT['fasta'].keys(),
            evalue = [blast_evalue]
        ) 
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast",
            "merge.annotation.blasn.tsv",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "merge_blastn.log"
        ),
    resources:
        cpus=1,
    conda:
        "../envs/biopython.yaml"
    threads: 1
    script: 
        "../scripts/merge_blastn.py"

