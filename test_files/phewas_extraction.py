import pandas as pd
from depricated.phewas_mp_functions import phewas_parallelization, combination_generation


def main():
    path = 'C:/Users/Sigve/Genome_Data/'
    input_file_names = [path + 'results/SNPs_non_coding.tsv',
                        path + 'results/SNPs_transcript_non_coding.tsv',
                        path + 'results/SNPs_missense.tsv',
                        path + 'results/SNPs_synonymous.tsv']

    snp_categories = ['non_coding', 'transcript_non_coding', 'missense', 'synonymous']

    snps = pd.DataFrame(columns=['variant_name', 'chrom', 'gene_name', 'gene_id'])
    category_count = 0

    for file in input_file_names:
        temp_df = pd.read_table(file)
        temp_df = temp_df[['variant_name', 'chrom', 'gene_name', 'gene_id']]
        temp_df['category'] = snp_categories[category_count]
        category_count += 1
        snps = pd.concat([snps, temp_df])

    snps.drop_duplicates(inplace=True)
    snps.reset_index(drop=True, inplace=True)

    phewas = pd.read_csv('C:/Users/Sigve/Genome_Data/SNP_data/phewas/phewas-catalog.csv')
    phewas['phewas code'] = phewas['phewas code'].apply(str)
    phewas = phewas[phewas['phewas code'].str.contains('290.\d*')]

    snps['phewas_codes'] = snps['variant_name'].apply(lambda x: set(phewas[phewas['snp'] == x]['phewas code'].tolist()))
    all_codes = set(phewas['phewas code'].tolist())

    df_list = []

    for code in all_codes:
        mask = snps['phewas_codes'].apply((lambda x: code in x))
        temp_df = snps[mask]

        temp_df['phewas_code'] = code

        if temp_df.empty:
            continue
        else:
            df_list.append(temp_df)

    combinations_df = phewas_parallelization(df_list, combination_generation, n_cores=12)
    combinations_df.to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/results/SNP_combinations_phewas_5_6.tsv', sep='\t')


if __name__ == '__main__':
    main()
