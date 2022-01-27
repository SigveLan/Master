import pandas as pd
import random


def main():
    """A script that generates a number of random combinations of up to 20 SNPs."""

    file_path = 'C:/Users/Sigve/Genome_Data/results/SNPs_non_synonymous.tsv'
    SNPs = pd.read_table(file_path, index_col=0)

    length = len(SNPs)
    result = {}
    for i in range(2000):
        sample = random.sample(range(0, length), random.randint(3, 9))

        # Gets the SNP IDs to a list.
        sample_SNPs = SNPs.filter(items=sample, axis=0).iloc[:, [0, 5]]

        sample_SNP_ids = ';'.join(set(sample_SNPs.iloc[:, 0].tolist()))
        # Gets the ensemble IDs to a list, and removes the version number.
        sample_gene_ids = ';'.join(set(sample_SNPs.iloc[:, 1].apply(lambda x: x.split('.')[0]).tolist()))
        result[i] = [sample_SNP_ids, sample_gene_ids]

    pd.DataFrame.from_dict(result, orient='index', columns=["SNP_ids", "gene_ids"]).to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/results/SNP_combinations.tsv', sep='\t')


if __name__ == '__main__':
    main()