import pandas as pd
import random


def main():

    file_path = 'C:/Users/Sigve/Genome_Data/results/SNPs_effect_mod.tsv'
    SNPs = pd.read_table(file_path, index_col=0)

    length = len(SNPs)
    result = {}
    for i in range(40):
        sample = random.sample(range(0, length), random.randint(1, 20))
        result[i] = ';'.join(SNPs.filter(items=sample, axis=0).iloc[:, 0].tolist())

    pd.Series(result, name='combinations').to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/results/SNP_combinations.tsv', sep='\t')


if __name__ == '__main__':
    main()