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
            "{sample}.selected.fasta"
        ),
    output:
        blast_out=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast",
            "virus",
            "{sample}",
            "{sample}.evalue_{evalue}.{database}.blastn.outfmt6.txt"
        ),
    params:
        database=lambda wildcards: os.path.join(
            DB_DICT["fasta"][wildcards.database]["path"],
            DB_DICT["fasta"][wildcards.database]["file"],
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "virus",
            "{sample}.evalue_{evalue}.{database}.blastn.outfmt6.log"
        ),    
    resources:
        cpus=2,
    conda:
        "../envs/blast.yaml"
    threads: 2
    shell:
        """
        if [ ! -e '{params.database}'.nsq ] && [ ! -e '{params.database}'.00.nsq ]; then
            makeblastdb -dbtype nucl -in '{params.database}' &> '{log}'
        fi

        blastn -query '{input.contig}' -db '{params.database}' -evalue {wildcards.evalue} \
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
    output:
        blast_out=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast",
            "human",
            "{sample}.nt.human.blastn.outfmt6.txt"
        ),
    params:
        database=blast_database,
        remote=blast_remote # -task blastn -remote
        threads=blast_threads_option # -num_threads {threads}
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "human",
            "{sample}.nt.human.blastn.outfmt6.log"
        ),    
    resources:
        cpus=2,
    conda:
        "../envs/blast.yaml"
    threads: blast_threads
    shell:
        """
        if [[ {params.database} != 'nt' ]] ; then
            if [ ! -e '{params.database}'.nsq ] && [ ! -e '{params.database}'.00.nsq ]; then
                makeblastdb -dbtype nucl -in '{params.database}' &> '{log}'
            fi
        fi

        blastn -query '{input.contig}' -db '{params.database}' -evalue 0.0001 \
                {params.remote} {params.threads} -taxids 9606\
               -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes  \
               -outfmt '6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle' \
               -min_raw_gapped_score 100 -perc_identity 90 -soft_masking true \
               -out '{output.blast_out}' -max_target_seqs 10 2>> '{log}'
        """