import pandas as pd
import random
import itertools


def main():
    """A script that generates a number of random combinations of a number of SNPs."""

    file_path = 'C:/Users/Sigve/Genome_Data/results/SNPs_for_FBA.tsv'
    SNPs = pd.read_table(file_path, index_col=0)

    length = len(SNPs)
    result = {}
    # Extensive is all combinations. If all combinations in a certain range is wanted, change the range below
    extensive = True

    if extensive:
        SNP_ids = SNPs['variant_name'].tolist()
        combs = []
        for L in range(2, length + 1):
            for subset in itertools.combinations(SNP_ids, L):
                combs.append(subset)
        gene_id_combs = []
        for subset in combs:
            gene_ids = []
            for i in subset:
                gene_ids.append(SNPs[SNPs['variant_name'] == i]['gene_id'].iloc[0])
            gene_id_combs.append(';'.join(gene_ids))

        combinations_df = pd.DataFrame({'snp_ids': combs, 'gene_ids': gene_id_combs})
        combinations_df['snp_ids'] = combinations_df['snp_ids'].apply(lambda x: ';'.join(x))
        combinations_df.to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/results/SNP_combinations.tsv', sep='\t')

    else:
        for i in range(20):
            sample = random.sample(range(0, length), random.randint(2, 4))

            # Gets the SNP IDs to a list.
            sample_SNPs = SNPs.filter(items=sample, axis=0).iloc[:, [0, 1]]

            sample_SNP_ids = ';'.join(set(sample_SNPs.iloc[:, 0].tolist()))
            # Gets the ensemble IDs to a list, and removes the version number.
            sample_gene_ids = ';'.join(set(sample_SNPs.iloc[:, 1].apply(lambda x: x.split('.')[0]).tolist()))
            result[i] = [sample_SNP_ids, sample_gene_ids]

        pd.DataFrame.from_dict(result, orient='index', columns=["snp_ids", "gene_ids"]).to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/results/SNP_combinations.tsv', sep='\t')


if __name__ == '__main__':
    main()