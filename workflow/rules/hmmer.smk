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