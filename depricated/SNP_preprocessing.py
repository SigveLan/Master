import pandas as pd
from src.assorted_functions import df_to_tsv

# Moved to assorted functions
path = 'C:/Users/Sigve/Genome_Data/'

SNP_file = path + 'SNP_data/dbSNP.tsv'
SNP_df = pd.read_table(SNP_file)

SNP_df = SNP_df[['#chr', 'pos', 'variation', 'snp_id']]
SNP_df = SNP_df[SNP_df['#chr'] != '#chr']
SNP_df.reset_index(inplace=True, drop=True)

SNP_df['variation'] = SNP_df['variation'].apply(lambda x: x.split('>'))
SNP_df[['ancestral', 'variation']] = pd.DataFrame(SNP_df['variation'].tolist())
SNP_df['variation'] = SNP_df['variation'].apply(lambda x: x.split(','))
SNP_df = SNP_df.explode('variation', ignore_index=True)
SNP_df['snp_id'] = SNP_df['snp_id'].apply(lambda x: 'rs' + str(x))
SNP_df['variation'] = SNP_df['ancestral'] + '/' + SNP_df['variation']
SNP_df.drop(columns='ancestral', inplace=True)
SNP_df = SNP_df[['snp_id', '#chr', 'pos', 'variation']]
SNP_df.rename(columns={'#chr': 'chr', 'variation': 'var'}, inplace=True)

SNP_df = SNP_df[SNP_df['var'].str.contains("^[ACTG]/[ACTG]$")]
SNP_df.sort_values(by=['chr'], inplace=True)
SNP_df.reset_index(inplace=True, drop=True)

df_to_tsv(SNP_df, path + 'SNP_data/dbSNP_edited.tsv')