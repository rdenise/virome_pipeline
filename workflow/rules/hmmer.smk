# Module containing all the hmmer related rules


##########################################################################
##########################################################################
# NOTES: 
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# Advantage of the ocncatenation, no nee to do it after by yourself
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule hmmpress:
    input:
        hmm=os.path.join(
            "{database_folder}",
            "{hmm_file}.hmm",
        ),
    output:
        h3i=os.path.join(
            "{database_folder}",
            "{hmm_file}.hmm.h3i",
        ),
        h3f=os.path.join(
            "{database_folder}",
            "{hmm_file}.hmm.h3f",
        ),
        h3m=os.path.join(
            "{database_folder}",
            "{hmm_file}.hmm.h3m",
        ),
        h3p=os.path.join(
            "{database_folder}",
            "{hmm_file}.hmm.h3p",
        ),
    resources:
        cpus=1,
    conda:
        "../envs/hmmer.yaml"
    threads: 1
    shell:
        """
        hmmpress '{input.hmm}'
        """


##########################################################################
##########################################################################

rule hmmsearch:
    input:
        proteins_fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "merge",
            "{sample}.faa",
        ),
        database=lambda wildcards: os.path.join(
            DB_DICT["hmm"][wildcards.database]["path"],
            DB_DICT["hmm"][wildcards.database]["file"],
        ),       
    output:
        tblout=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "hmmer",
            "{sample}",
            "{sample}.{database}.tblout.txt",
        ),
        domtblout=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "hmmer",
            "{sample}",
            "{sample}.{database}.domtblout.txt",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "hmmer",
            "{sample}.{database}.hmmsearch.log"
        ),  
    resources:
        cpus=5,
    conda:
        "../envs/hmmer.yaml"          
    threads: 5     
    shell:
        """
        hmmsearch --domtblout '{output.domtblout}' --tblout '{output.tblout}'\
         -E 1.0 --cpu {threads} '{input.database}' '{input.proteins_fasta}' &> '{log}'
        """


##########################################################################
##########################################################################


rule concat_results_domtblout :
    input : 
        expand(os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "hmmer",
            "{sample}", 
            "{sample}.{database}.domtblout.txt"
            ),
        sample=CONTIGS_DICT.keys(),
        database = DB_DICT['hmm'].keys(),
        )
    output :
        domtblout=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "hmmer",
            "merge.domtblout.txt"
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "hmmer",
            "concat_results_domtblout.log"
        ),  
    resources:
        cpus=1,         
    threads : 1    
    shell :
        """
        echo 'contig_id\taccession_target\ttlen\tgene\taccession_query\tqlen\tE_value_full\tscore_full\tbias_full\tdom_num\ttotal_dom_hit\tc_Evalue\ti_Evalue\tdom_score\tbias_dom\tq_start\tq_stop\tali_start\tali_stop\tenv_start\tenv_stop\tacc' > '{output}'

        for domtblout in {input}
        do
            if grep -q ^[^#] "$domtblout"
            then 
                grep ^[^#] "$domtblout" | sed -E $'s/ +/\t/g' | cut -d $'\t' -f 1-22 >> '{output}' 2> '{log}'
            fi
        done
        """ 

##########################################################################
##########################################################################


rule concat_results_tblout :
    input : 
        expand(os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "hmmer",
            "{sample}", 
            "{sample}.{database}.tblout.txt"
            ),
        sample=CONTIGS_DICT.keys(),
        database = DB_DICT['hmm'].keys(),
        )
    output :
        tblout=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "hmmer",
            "merge.tblout.txt"
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "hmmer",
            "concat_results_tblout.log"
        ),  
    resources:
        cpus=1,
    threads : 1    
    shell :
        """
        echo 'contig_ig\taccession_target\tgene\taccession_query\tE_value_full\tscore_full\tbias_full\tE-value_dom\tscore_dom\tbias_dom\texp\treg\tclu\tov\tenv\tdom\trep\tinc' > '{output}'

        for tblout in {input}
        do
            if grep -q ^[^#] "$tblout"
            then
                grep ^[^#] "$tblout" | sed -E $'s/ +/\t/g' | cut -d $'\t' -f 1-18 >> '{output}' 2> '{log}'
            fi
        done
        """    


##########################################################################
##########################################################################
