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


rule virsorter_postprocess_step1:
    input:
        final_score=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "vs2-step1",
            "final-viral-score.tsv",
        ),
        contamination=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "{sample}",
            "contamination.tsv",
        ),
    output:
        ids_keep1=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter_positive.keep1.ids",
        ),
        ids_keep2=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter_positive.keep2.ids",
        ),
        manual_check=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter2check.ids",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "virsorter",
            "{sample}.virsorter_postprocess_step1.log",
        ),
    resources:
        cpus=1,
    conda:
        "../envs/pandas_plots.yaml"
    threads: 1
    script: 
        "../scripts/virsorter_postprocess_step1.py"


##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule virsorter_postprocess_step2:
    input:
        ids_keep2=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter_positive.keep2.ids",
        ),
        tsv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "annotations.tsv"
        ),
        suspicous_gene="config/suspicious-gene.list",
    output:
        suspicious=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter_positive.keep2.suspicious.ids",
        ),
        checked=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter_positive.keep2.checked.ids",
        ),
        manual_check=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter2check.ids",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "virsorter",
            "{sample}.virsorter_postprocess_step2.log",
        ),
    resources:
        cpus=1,
    conda:
        "../envs/pandas_plots.yaml"
    threads: 1
    script: 
        "../scripts/virsorter_postprocess_step2.py"


##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule combine_virsorter_virfinder:
    input:
        ids_virsorter_keep2=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter_positive.keep2.checked.ids",
        ),
        ids_virsorter_keep1=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter_positive.keep1.ids",
        ),
        ids_virfinder=lambda wildcards: os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "deepvirfinder",
            wildcards.sample,
            f"{wildcards.sample}.deepvirfinder_positive.gt{cutoff_deepvirfinder}bp.ids",
        ),  
        contigs=lambda wildcards: os.path.join(
            CONTIGS_FOLDER,
            CONTIGS_DICT[wildcards.sample]["file"],
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
                "virus",
                "{sample}.evalue_{evalue:.0e}.{database}.blastn.outfmt6.txt"
            ),
            sample = CONTIGS_DICT.keys(),
            database = DB_DICT['fasta'].keys(),
            evalue = [blast_evalue],
        ) 
    output:
        tsv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast",
            "virus",
            "merge.eval_{evalue}.cov_{coverage}.annotation.blasn.tsv",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "virus",
            "merge_blastn.eval_{evalue}.cov_{coverage}.log"
        ),
    resources:
        cpus=1,
    conda:
        "../envs/biopython.yaml"
    threads: 1
    script: 
        "../scripts/merge_blastn.py"


##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule merge_blastn_human:
    input:
        all_out=expand(
            os.path.join(
                OUTPUT_FOLDER,
                "processing_files",
                "blast",
                "human",
                "{sample}.evalue_{evalue:.0e}.{database}.human.blastn.outfmt6.txt"
            ),
            sample = CONTIGS_DICT.keys(),
            database = DB_DICT['human'].keys(),
            evalue = [blast_evalue],
        ) 
    output:
        tsv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast",
            "human",
            "merge.eval_{evalue}.cov_{coverage}.human.annotation.blasn.tsv",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "human",
            "merge_blastn.eval_{evalue}.cov_{coverage}.log"
        ),
    resources:
        cpus=1,
    conda:
        "../envs/biopython.yaml"
    threads: 1
    script: 
        "../scripts/merge_blastn.py"

##########################################################################
##########################################################################

rule postprocess_hmmsearch :
    input :
        tblout=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "hmmer",
            "merge.tblout.txt"
            ),
        domtblout=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "hmmer",
            "merge.domtblout.txt"
            ),
    output :
        significant_hit = os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "hmmer",
            f"significant_hit.full_{hmm_evalue_full:.0e}.dom_{hmm_evalue_dom}.domtblout.txt"
            ),    
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "postprocess_hmmsearch.log"
        ),
    resources:
        cpus=1,
    conda:
        "../envs/pandas_plots.yaml"
    threads: 1
    script: 
        "../scripts/postprocess_hmmsearch.py"


##########################################################################
##########################################################################
