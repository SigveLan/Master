import pandas as pd
import re
# A simple script to sort SNPs by position, and also filter out duplicates. Not really necessary


def main():

    SNPs_file = 'C:/Users/Sigve/Genome_Data/SNP_data/biomart_chrom_all_test_orig.tsv'
    SNP_df = pd.read_table(SNPs_file)
    SNP_df.drop_duplicates(subset=['Variant name'], inplace=True)
    SNP_df.sort_values(by=['Chromosome/scaffold position start (bp)'], inplace=True)
    SNP_df.set_index(['Variant name'], drop=True, inplace=True)
    SNP_df['Variant alleles'] = SNP_df['Variant alleles'].apply(str)
    SNP_df = SNP_df[SNP_df['Variant alleles'].str.contains("^[ACTG]/[ACTG]$")]

    SNP_df['Chromosome/scaffold position start (bp)'] = SNP_df['Chromosome/scaffold position start (bp)'].apply(int)
    SNP_df['Chromosome/scaffold position end (bp)'] = SNP_df['Chromosome/scaffold position end (bp)'].apply(int)
    SNP_df['Strand'] = SNP_df['Strand'].apply(int)



    SNP_df.to_csv('C:/Users/Sigve/Genome_Data/SNP_data/biomart_chrom_all_test.tsv', sep='\t')


if __name__ == '__main__':
    main()
