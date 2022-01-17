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
            OUTPUT_FOLDER,
            "databases",
            "viral_contigs",
            "{sample}.selected.fasta",
        ),
        database=lambda wildcards: os.path.join(
            DB_DICT["fasta"][wildcards.database]["path"],
            DB_DICT["fasta"][wildcards.database]["file"],
        ),
    output:
        blast_out=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast",
            "{sample}.evalue_{evalue}.{database}.blastn.outfmt6.txt"
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "{sample}.evalue_{evalue}.{database}.blastn.outfmt6.log"
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
        if [ ! -e '{input.database}'.ndb ]; then
            makeblastdb -dbtype nucl -in '{input.database}' &> '{log}'
        fi

        blastn -query '{input.contig}' -db '{input.database}' -evalue {wildcards.evalue} \
               -outfmt '6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle' \
               -out '{output.blast_out}' -num_threads {threads} -num_alignments 25000 2>> '{log}'
        """


##########################################################################
##########################################################################


rule blastn_human:
    input:
        contig=lambda wildcards: os.path.join(
            CONTIGS_FOLDER,
            CONTIGS_DICT[wildcards.sample]["file"],
        ),
        database=lambda wildcards: os.path.join(
            DB_DICT["human"]["nt"]["path"],
            DB_DICT["human"]["nt"]["file"],
        ),
    output:
        blast_out=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast",
            "{sample}.evalue_{evalue}.nt.blastn.outfmt6.txt"
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "{sample}.evalue_{evalue}.nt.blastn.outfmt6.log"
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
        if [ ! -e '{input.database}'.ndb ]; then
            makeblastdb -dbtype nucl -in '{input.database}' &> '{log}'
        fi

        blastn -query '{input.contig}' -db '{input.database}' -evalue {wildcards.evalue} \
               -outfmt '6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle' \
               -out '{output.blast_out}' -num_threads {threads} -taxids 9606 -num_alignments 25000 2>> '{log}'
        """