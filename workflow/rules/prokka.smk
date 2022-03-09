# Module containing all the ncbi-blast related rules

##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# Advantage of the concatenation, no nee to do it after by yourself
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule prokka:
    input:
        contig=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "viral_contigs",
            "{sample}.selected.fasta",
        ),
        database_blast=prokka_protein_db,
        h3i=prokka_hmm_db + ".h3i",       
    output:
        output_dir=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "prokka",
            "{sample}",
            "prokka_name",
            "{sample}.prokka.pvogs.crass.tbl",
        ),
    params:
        output_dir=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "prokka",
            "{sample}",
            "prokka_name",
        ),
        prefix="{sample}.prokka.pvogs.crass",
        gcode=11,
        hmm=prokka_hmm_db,
        kingdom=prokka_kingdom,      
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "prokka",
            "{sample}.log"
        ),
    resources:
        cpus=5,
    conda:
        "../envs/prokka.yaml"
    threads: 5
    shell:
        """
        prokka --outdir '{params.output_dir}' --prefix '{params.prefix}' \
        --gcode '{params.gcode}' --hmms '{params.hmm}'  \
        --proteins '{input.database_blast}' \
        --compliant --partialgenes --cpus '{threads}' \
        --kingdom '{params.kingdom}' '{input.contig}' --force &> '{log}'
        """


##########################################################################
##########################################################################

rule prokka_rename :
    input :
        tbl=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "prokka",
            "{sample}",
            "prokka_name",
            "{sample}.prokka.pvogs.crass.tbl",
        ),
        fasta_contig=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "viral_contigs",
            "{sample}.selected.fasta",
        ),
    output :
        faa=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "prokka",
            "{sample}",
            "contig_renamed",
            "{sample}.prokka.pvogs.crass.faa",
        ),
        transTbl=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "prokka",
            "{sample}",
            "{sample}.prokka.pvogs.crass.prokka_locustag.tsv",
        ),        
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "prokka",
            "{sample}.rename.log"
        ),
    resources:
        cpus=1,
    conda:
        "../envs/biopython.yaml"
    threads: 1
    script: 
        "../scripts/rename_prokka.py"


##########################################################################
##########################################################################
