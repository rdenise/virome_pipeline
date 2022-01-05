# Module containing all the ncbi-blast related rules


##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# Advantage of the ocncatenation, no nee to do it after by yourself
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule blastn:
    input:
        contig=os.path.join(
            CONTIG_FOLDER, 
            "??_solo_contig_or_merge_contig??.fasta"
        ),
        database_blast=os.path.join(
            DATABASE_FOLDER, 
            "{database}", 
            "viral_proteins_plus_crass.fasta"
        ),
        h3i=os.path.join(
            DATABASE_FOLDER,
            "all_vogs.hmm.h3i"
        ),        
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "prokka",
            "{contig}",
            "{contig}.proteins.fa"
        ),
    params:
        output_dir=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "prokka",
            "{contig}",
        ),
        prefix="{contig}.prokka.pvogs.crass",
        gcode=11,
        hmm=os.path.join(
            DATABASE_FOLDER,
            "all_vogs.hmm"
        ),         
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "prokka",
            "contig_{contig}.log"
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
        --locus-tag PALEO --compliant
        '{input.contig}' &> '{log}'
        """


##########################################################################
##########################################################################
