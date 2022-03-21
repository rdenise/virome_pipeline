import pandas as pd
import sys

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################
# From the method for virsorter 2:
# Keep1: viral_gene >0
# Keep2: viral_gene =0 AND (host_genes =0 OR score >=0.95 OR hallmark >2)
# Manual check: (NOT in Keep1 OR Keep2) AND viral_gene =0 AND host_genes =1 AND length >=10kb
# Discard: the rest

# To look at the viral_gene, host_genes, score, and hallmark of sequences you can merge "vs2-pass1/final-viral-score.tsv" and "checkv/contamination.tsv", and filter in a spreadsheet.

# Get the informations from virsorter
virsorter_step1 = snakemake.input.final_score
virsorter_df = pd.read_table(virsorter_step1)
virsorter_df = virsorter_df.rename(columns={"seqname":"contig_id"})

# Get the information from CheckV
checkv = snakemake.input.contamination
checkv_df = pd.read_table(checkv)

# Need to check what is the output of checkv
merge_df = virsorter_df.merge(checkv_df, how='outer', on='contig_id')

# First filter: viral_gene >0
keep1_contig_ids = merge_df[merge_df.viral_genes > 0].contig_id.unique().tolist()


# Write in a file
with open(snakemake.output.ids_keep1, "w") as w_file:
    w_file.write("contig_id\n")

    for name in keep1_contig_ids:
        w_file.write(f"{name}\n")

# Removing Keep1
merge_tmp = merge_df[~(merge_df.contig_id.isin(keep1_contig_ids))]


# Second filter: viral_gene =0 AND (host_genes =0 OR score >=0.95 OR hallmark >2)
keep2_contig_ids =  merge_tmp[(merge_tmp.viral_genes == 0) & \
                          ((merge_tmp.host_genes == 0) | \
                           (merge_tmp.max_score >= 0.95) | \
                           (merge_tmp.hallmark > 2))].contig_id.unique().tolist()

# Write in a file
with open(snakemake.output.ids_keep2, "w") as w_file:
    w_file.write("contig_id\n")

    for name in keep2_contig_ids:
        w_file.write(f"{name}\n")

# Removing Keep2
merge_tmp = merge_tmp[~(merge_tmp.contig_id.isin(keep2_contig_ids))]

# Manualcheck: (NOT in Keep1 OR Keep2) AND viral_gene =0 AND host_genes =1 AND length >=10kb
manualcheck_contig_ids = merge_tmp[(merge_tmp.viral_genes == 0) & \
                                (merge_tmp.host_genes == 1) & \
                                (merge_tmp.length >= 10000)].contig_id.tolist()


# Write the big dataframe to combine with DRAMv annotation later
merge_tmp.loc[(merge_tmp.contig_id.isin(manualcheck_contig_ids)),:].to_csv(snakemake.output.manual_check_tsv, sep="\t", index=False)

# Write the ids in a file
merge_tmp.loc[(merge_tmp.contig_id.isin(manualcheck_contig_ids)),'contig_id'].to_csv(snakemake.output.manual_check_ids, sep="\t", index=False)

# Write discarded contig
merge_tmp.loc[~(merge_tmp.contig_id.isin(manualcheck_contig_ids)),:].to_csv(snakemake.output.discarded_tsv, sep="\t", index=False)

# Write the ids in a file
merge_tmp.loc[~(merge_tmp.contig_id.isin(manualcheck_contig_ids)),'contig_id'].to_csv(snakemake.output.discarded_ids, sep="\t", index=False)


###########################################################
###########################################################

