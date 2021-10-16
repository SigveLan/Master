import pandas as pd

"""A script to add gene number from the CBM model to the SNP results"""


def main():

    model_file = 'C:/Users/Sigve/Genome_Data/exon_model_data/recon-store-genes.tsv'
    model_df = pd.read_table(model_file)
    model_df = model_df[['gene_number', 'ensembl_gene']]

    model_df['gene_number'] = model_df['gene_number'].astype(int)
    model_df.drop_duplicates(subset=['gene_number'], inplace=True)

    filtered_SNPs_file = 'C:/Users/Sigve/Genome_Data/results/SNPs_effect.tsv'
    SNP_results = pd.read_table(filtered_SNPs_file, index_col=0)

    SNP_results.set_index(SNP_results['gene_id'].apply(lambda x: str(x).split('.')[0]), inplace=True)

    joined_df = SNP_results.join(model_df.set_index('ensembl_gene'))

    joined_df.dropna(subset=['gene_number'], inplace=True)
    joined_df.drop_duplicates(subset=['amino_acid_change'], inplace=True)
    joined_df.reset_index(drop=True, inplace=True)
    joined_df['gene_number'] = joined_df['gene_number'].astype(int)

    joined_df.to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/results/SNPs_effect_mod.tsv', sep='\t')


if __name__ == '__main__':
    main()