import sys
import os
import pandas as pd

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

domtblout = snakemake.input.domtblout

columns = ['contig_id', 'accession_target', 'tlen', 'gene', 
           'accession_query', 'qlen', 'E_value_full', 'score_full', 
           'bias_full', 'dom_num', 'total_dom_hit', 'c_Evalue', 
           'i_Evalue', 'dom_score', 'bias_dom', 'q_start', 
           'q_stop', 'ali_start', 'ali_stop', 'env_start', 'env_stop', 'acc']

dict_columns = {column:columns.index(column) for column in columns}

with open(snakemake.output.significant_hit, 'w') as w_file:
    with open(domtblout) as r_file:
        rstrip_line = r_file.readline().rstrip()
        w_file.write(f"{rstrip_line}\tdatabase\n")

        for line in r_file:
            split_line = line.split()

            # profile_coverage = split_line[dict_columns['q_stop']] - split_line[dict_columns['q_start']] / split_line[dict_columns['qlen']]
            # sequence_coverage = split_line[dict_columns['env_stop']] - split_line[dict_columns['env_start']] / split_line[dict_columns['tlen']]
        
            # coverage = min(profile_coverage, sequence_coverage)

            if split_line[dict_columns['i_Evalue']] < snakemake.wildcards.hmm_evalue_dom and split_line[dict_columns['E_value_full']] < snakemake.wildcards.hmm_evalue_full:
                rstrip_line = line.rstrip()
                database = "pVOGs" if split_line[dict_columns['gene']].startswith('VOG') else PHROGs
                w_file.write(f"{rstrip_line}\t{database}\n")

###########################################################
###########################################################

