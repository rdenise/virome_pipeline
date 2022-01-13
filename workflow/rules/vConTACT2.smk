##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule vcontact2_preprocess:
    input:
        proteins_fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "prokka",
            "{contig}",
            "{contig}.prokka.pvogs.crass.faa",
        ),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "prodigal"
            "{contig}_translations_vir_contigs1000.rename.faa"
        ),
        csv=os.path.join(
            OUTPUT_FOLDER,
            "prodigal"
            "{contig}_proteins.csv"
        ),        
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "vcontact2",
            "contig_{contig}.vcontact2_preprocess.log"
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
            "prodigal"
            "{contig}_translations_vir_contigs1000.rename.faa"
        ),
        protein_csv=os.path.join(
            OUTPUT_FOLDER,
            "prodigal"
            "{contig}_proteins.csv"
        ),
    output:
        csv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "vcontact2",
            "contig_{contig}",
            "genome_by_genome_overview.csv",
        ),
    params:
        vcontact2_db="ProkaryoticViralRefSeq85-Merged",
        rel_mode="Diamond",
        pcs_mode="MCL",
        vcs_mode="ClusterONE",
        c1_bin="cluster_one-1.0.jar",
        output_dir=directory(os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "vcontact2",
            "contig_{contig}"
        )),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "vcontact2",
            "contig_{contig}.vcontact2.log"
        ),
    resources:
        cpus=5,
    conda:
        "../envs/vcontact2.yaml"
    threads: 5
    shell:
        """
        vcontact2 --raw-proteins '{input.fasta}' --rel-mode {params.rel_mode} \
        --proteins-fp '{input.protein_csv}' --db {params.vcontact2_db} \
        --pcs-mode {params.pcs_mode} --vcs-mode {params.vcs_mode} --c1-bin {params.c1_bin} \
        --output-dir {params.output_dir} -t {threads}
        """


##########################################################################
##########################################################################
