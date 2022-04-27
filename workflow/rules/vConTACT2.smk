##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being
# Important:: Need to maybe separate the contig inside the BIG contig file because each could be 
# a bacteria or a virus by itself

rule vcontact2_preprocess:
    input:
        proteins_fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "merge",
            "{sample}.faa",
        ),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "merge",
            "{sample}.rename.faa",
        ),
        csv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "merge",
            "{sample}.proteins.csv"
        ),
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "merge",
            "{sample}.rename.low.faa",
        ),
    params:
        fasta_contig=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "viral_contigs",
            "{sample}.selected.fasta",
        ),     
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "vcontact2",
            "contig_{sample}.vcontact2_preprocess.log"
        ),
    resources:
        cpus=1,
    conda:
        "../envs/biopython.yaml"
    threads: 1
    script:
        "../scripts/protein_csv.py"


##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule vcontact2:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "merge",
            "{sample}.rename.faa",
        ),
        protein_csv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "merge",
            "{sample}.proteins.csv"
        ),
    output:
        csv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "vcontact2",
            "{sample}",
            "genome_by_genome_overview.csv",
        ),
    params:
        vcontact2_db="ProkaryoticViralRefSeq207-Merged",
        rel_mode="Diamond",
        pcs_mode="MCL",
        vcs_mode="ClusterONE",
        c1_bin="cluster_one-v1.0.jar",
        output_dir=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "vcontact2",
            "{sample}",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "vcontact2",
            "{sample}.vcontact2.log"
        ),
    resources:
        cpus=5,
    conda:
        "../envs/vcontact2.yaml"
    threads: 5
    shell:
        """
        findClusterONE=( $(find . -name '{params.c1_bin}') )

        vcontact2 --raw-proteins '{input.fasta}' --rel-mode {params.rel_mode} \
        --proteins-fp '{input.protein_csv}' --db {params.vcontact2_db} \
        --pcs-mode {params.pcs_mode} --vcs-mode {params.vcs_mode} \
        --c1-bin "${{findClusterONE[0]}}" \
        --output-dir {params.output_dir} -t {threads} &> '{log}'
        """


##########################################################################
##########################################################################
