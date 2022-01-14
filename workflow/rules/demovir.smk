##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being
# 3. https://github.com/feargalr/Demovir need to be recoded to be better used in conda or in the pipeline
# Or can be "recoded" here because it doesn't seems difficult from the script I saw


rule demovir:
    input:
        contig=os.path.join(
            CONTIG_FOLDER, 
            "AA.fasta"
        ),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "demovir",
            "{sample}",
            "all.fa_circular.fna",
        ),        
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "demovir",
            "{sample}.demovir.log"
        ),
    resources:
        cpus=1,
    conda:
        "../envs/virsorter.yaml"
    threads: 1
    shell:
        """
        demovir or I recode the thing because needed
        """


##########################################################################
##########################################################################
