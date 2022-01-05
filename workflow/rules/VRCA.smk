##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being
# 3. VRCA (before named VICA) need to be install on a conda environment first because it is only a script
# Also the information could be extract from VirSorter probably
# Or can be "recoded" here because it doesn't seems difficult from the script I saw
# https://github.com/alexcritschristoph/VRCA


rule vrca:
    input:
        contig=os.path.join(
            CONTIG_FOLDER, "??_solo_contig_or_merge_contig??.fasta"
        ),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "VRCA",
            "contig_{contig}",
            "all.fa_circular.fna",
        ),        
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "VRCA",
            "contig_{contig}.VRCA.log"
        ),
    resources:
        cpus=1,
    conda:
        "../envs/virsorter.yaml"
    threads: 1
    shell:
        """
        python find_circular.py -i all.fa
        """


##########################################################################
##########################################################################
