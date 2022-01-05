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
            CONTIG_FOLDER, "??_solo_contig_or_merge_contig??.fasta"
        ),
        database=os.path.join(
            DATABASE_FOLDER, "{database}", "??.fasta"
        ),
    output:
        blast_out=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast"
            "contig_{contig}.evalue_1e-10.{database}.blastn.outfmt6.txt"
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "contig_{contig}.evalue_1e-10.{database}.blastn.outfmt6.log"
        ),
    resources:
        cpus=5,
    conda:
        "../envs/blast.yaml"
    envmodules:
        "ncbi_blast/2.10.1",
    threads: 5
    shell:
        """
        makeblastdb -dbtype prot -in '{input.database}' &> '{log}'

        blastn -query '{input.contig}' -db '{input.database}' -evalue 1e-10 \
               -outfmt '6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle' \
               -out '{output.blast_out}' -num_threads {threads} -num_alignments 25000 2>> '{log}'
        """


##########################################################################
##########################################################################
