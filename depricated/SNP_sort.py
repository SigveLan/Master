import pandas as pd

# A simple script to sort SNPs by position. Not really necessary

SNPs_file = 'C:/Users/sigve/Documents/Genome_Data/SNP_data/biomart_chrom_1_test_1.tsv'

SNP_df = pd.read_table(SNPs_file, index_col=0)
SNP_df.sort_values(by=['Chromosome/scaffold position start (bp)'], inplace=True)
SNP_df.to_csv('C:/Users/sigve/Documents/Genome_Data/SNP_data/biomart_chrom_1_test_1.tsv', sep='\t')