import pandas as pd
from src.assorted_functions import df_to_tsv

# Filters out any SNPs that will lead to multiple FBA runs of the same deletions.

path = 'C:/Users/Sigve/Genome_Data/'

SNP_file = path + 'results/SNPs_missense.tsv'
SNP_df = pd.read_table(SNP_file)

# Keep first duplicate value
SNP_df = SNP_df.drop_duplicates(subset=['variant_name'])
SNP_df.reset_index(inplace=True, drop=True)
SNP_df['gene_id'] = SNP_df['gene_id'].apply(lambda x: x.split('.')[0])

df_to_tsv(SNP_df[['variant_name', 'gene_id', 'score']], path + 'results/SNPs_for_FBA.tsv')
