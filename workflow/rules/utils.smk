# Module with utilitary module as untar, ...


##########################################################################
##########################################################################

rule uncompress:
    input:
        contig=lambda wildcards: os.path.join(
            CONTIGS_FOLDER,
            CONTIGS_DICT[wildcards.contig]["file"] + \
            CONTIGS_DICT[wildcards.contig]["ext_compress"],
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
            "{contig}.{ext}.log"
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
            CONTIGS_DICT[wildcards.contig]["path"],
            CONTIGS_DICT[wildcards.contig]["file"],
        ),
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "clean",
            "{contig}.log"
        ),    
    threads: 1
    shell:
        """
        rm "{input.contig}"
        """


##########################################################################
##########################################################################
