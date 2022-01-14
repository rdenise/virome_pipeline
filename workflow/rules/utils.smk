# Module with utilitary module as untar, ...


##########################################################################
##########################################################################

rule uncompress:
    input:
        contig=lambda wildcards: os.path.join(
            CONTIGS_FOLDER,
            CONTIGS_DICT[wildcards.sample]["file"] + \
            CONTIGS_DICT[wildcards.sample]["ext_compress"],
        ),
    output:
        contig = os.path.join(
            CONTIGS_FOLDER,
            "{contig,[^.]+}.{ext}",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "uncompress",
            "{sample}.{ext}.log"
        ),    
    threads: 1
    shell:
        """
        if [[ "{input.contig}" == *.tar.gz ]]; then
            tar -xzvf "{input.contig}"
        else 
            gzip -dk "{input.contig}"
        fi
        """


##########################################################################
##########################################################################

rule clear_uncompress:
    input:
        contig=lambda wildcards: os.path.join(
            CONTIGS_DICT[wildcards.sample]["path"],
            CONTIGS_DICT[wildcards.sample]["file"],
        ),
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "clean",
            "{sample}.log"
        ),    
    threads: 1
    shell:
        """
        rm "{input.contig}"
        """


##########################################################################
##########################################################################
