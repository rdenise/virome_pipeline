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
            "{sample}.{cutoff}.deepvirfinder_postprocess.log",
        ),
    resources:
        cpus=1,
    threads: 1
    script:
        "../scripts/deepvirfinder_postprocess.py"


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
        manual_check_tsv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter.need_manual_check.tsv",
        ),
        discarded_tsv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter.discarded.tsv",
        ),
        manual_check_ids=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter.need_manual_check.ids",
        ),
        discarded_ids=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter.discarded.ids",
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
        annotations=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "annotations.tsv",
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


rule extract_deepvirfinder:
    input:
        ids_virfinder=lambda wildcards: os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "deepvirfinder",
            wildcards.sample,
            f"{wildcards.sample}.deepvirfinder_positive.gt{cutoff_deepvirfinder}bp.ids",
        ),
        contig=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "contigs",
            "human_filtered",
            "{sample}.filtered.sorted.fasta",
        ),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "deepvirfinder",
            "{sample}",
            "{sample}.extract.fasta",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "deepvirfinder",
            "{sample}.extract_deepvirfinder.log",
        ),
    resources:
        cpus=1,
    conda:
        "../envs/pandas_plots.yaml"
    threads: 1
    script:
        "../scripts/extract_deepvirfinder.py"


##########################################################################
##########################################################################
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule combine_virsorter_virfinder:
    input:
        ids_virsorter_keep2_checked=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter_positive.keep2.checked.ids",
        ),
        ids_virsorter_keep2_suspicious=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter_positive.keep2.suspicious.ids",
        ),
        ids_virsorter_manual_check=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter.need_manual_check.ids",
        ),
        ids_virsorter_discarded=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "virsorter.discarded.ids",
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
        contigs_virsorter=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "final-viral-combined.fa",
        ),
        contigs_deepvirfinder=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "deepvirfinder",
            "{sample}",
            "combined.fa",
        ),        
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "viral_contigs",
            "{sample}.selected.fasta",
        ),
        tsv=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "viral_contigs",
            "{sample}.selected.tsv",
        ),
        translation_table=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "viral_contigs",
            "{sample}.translation_table_contig.tsv",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "postprocess_detection",
            "{sample}.combine_virsorter_virfinder.log",
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


rule prepare_missing_annotation:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "viral_contigs",
            "{sample}.selected.fasta",
        ),
        tsv=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "viral_contigs",
            "{sample}.selected.tsv",
        ),
        translation_table=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "viral_contigs",
            "{sample}.translation_table_contig.tsv",
        ),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "viral_contigs",
            "{sample}.missing_annotation.fasta",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "postprocess_detection",
            "{sample}.prepare_missing_annotation.log",
        ),
    resources:
        cpus=1,
    conda:
        "../envs/biopython.yaml"
    threads: 1
    script:
        "../scripts/prepare_missing_annotation.py"


##########################################################################
##########################################################################


rule merge_annotations:
    input:
        tsv_missing=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "missing_annotation",
            "annotations.tsv",
        ),
        faa_missing=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "missing_annotation",
            "{sample}.faa",
        ),
        tsv_virsorter=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "annotations.tsv",
        ),
        faa_virsorter=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "{sample}.faa",
        ),
        translation_table=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "viral_contigs",
            "{sample}.translation_table_contig.tsv",
        ),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "merge",
            "{sample}.faa",
        ),
        tsv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "merge",
            "annotations.tsv",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "postprocess_detection",
            "{sample}.merge_annotations.log",
        ),
    resources:
        cpus=1,
    conda:
        "../envs/biopython.yaml"
    threads: 1
    script:
        "../scripts/merge_annotations.py"


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
                "{sample}",
                "{sample}.evalue_{evalue:.0e}.{database}.blastn.outfmt6.txt",
            ),
            sample=CONTIGS_DICT.keys(),
            database=DB_DICT["fasta"].keys(),
            evalue=[blast_evalue],
        ),
    output:
        tsv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast",
            "virus",
            "merge.eval_{evalue}.cov_{coverage}.pident_{pident}.annotation.blasn.tsv",
        ),
    params:
        contigs=expand(
            os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "viral_contigs",
                "{sample}.selected.fasta",
            ),
            sample=CONTIGS_DICT.keys(),
        ),
        dict_databases=DB_DICT["fasta"],
        minimum_length=config["default_blast_option"]["length_min"],
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "virus",
            "merge_blastn.eval_{evalue}.cov_{coverage}.pid_{pident}.log",
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
        tsv=lambda wildcards: os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast",
            "human",
            f"{wildcards.sample}.nt.human.blastn.outfmt6.txt",
        ),
        fasta=lambda wildcards: os.path.join(
            CONTIGS_FOLDER,
            CONTIGS_DICT[wildcards.sample]["file"],
        ),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "contigs",
            "human_filtered",
            "{sample}.filtered.sorted.fasta",
        ),
    params:
        tsv=lambda wildcards: os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast",
            "human",
            f"{wildcards.sample}.cov_0.6.human.annotation.blasn.tsv",
        ),
        minimum_length=config["default_blast_option"]["length_min"],
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "human",
            "{sample}_blastn.filtered.human.log",
        ),
    resources:
        cpus=1,
    conda:
        "../envs/biopython.yaml"
    threads: 1
    script:
        "../scripts/remove_human.py"


##########################################################################
##########################################################################


rule postprocess_hmmsearch:
    input:
        tblout=os.path.join(
            OUTPUT_FOLDER, "processing_files", "hmmer", "merge.tblout.txt"
        ),
        domtblout=os.path.join(
            OUTPUT_FOLDER, "processing_files", "hmmer", "merge.domtblout.txt"
        ),
    output:
        significant_hit=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "hmmer",
            "significant_hit.full_{hmm_evalue_full}.dom_{hmm_evalue_dom}.domtblout.txt",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "postprocess_hmmsearch.full_{hmm_evalue_full}.dom_{hmm_evalue_dom}.log",
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
