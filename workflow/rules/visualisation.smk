##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule annotate_viral_taxonomy:
    input:
        merge_blast=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast",
            "virus",
            f"merge.eval_{blast_evalue:.0e}.cov_{blast_coverage}.pident_{blast_pident}.annotation.blasn.tsv"
        ),
    output:
        tsv = os.path.join(
            OUTPUT_FOLDER,
            "results",
            "taxonomic_annotation_contigs.tsv"
        ),
    params:
        metadata=config["annotation_phages"],
        viral_tsv=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "viral_contigs",
        ),
        virsorter_tsv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "annotation",
            "annotate_viral_taxonomy.log"
        ),
    resources:
        cpus=1,
    conda:
        "../envs/biopython.yaml"
    threads: 1
    script: 
        "../scripts/contig_annotation_blast.py"


##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule annotate_viral_protein:
    input:
        merge_hmmout=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast",
            "virus",
            f"merge.eval_{blast_evalue:.0e}.cov_{blast_coverage}.pident_{blast_pident}.annotation.blasn.tsv"
        ),
    output:
        tsv = os.path.join(
            OUTPUT_FOLDER,
            "results",
            "taxonomic_annotation_contigs.tsv"
        ),
    params:
        metadata=config["annotation_phages"],
        viral_tsv=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "viral_contigs",
        ),
        virsorter_tsv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "annotation",
            "annotate_viral_taxonomy.log"
        ),
    resources:
        cpus=1,
    conda:
        "../envs/biopython.yaml"
    threads: 1
    script: 
        "../scripts/contig_annotation_blast.py"

