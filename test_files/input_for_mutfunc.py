import pandas as pd


def main():

    path = 'C:/Users/Sigve/Genome_Data/'
    input_file = path + 'results/SNPs_missense.tsv'
    output_file = path + 'results/SNPs_mutfunc.tsv'

    df = pd.read_table(input_file, index_col=0)
    df['amino_acid_pos'] = df['amino_acid_pos'].astype(int).astype(str)

    df = df[df['amino_acid_change'].map(lambda x: x[2] != '_')]
    df[['gene_name', 'amino_acid_pos', 'amino_acid_change']].to_csv(output_file, sep='\t', index=False)

    #df[['chrom', 'chrom_pos', 'variant_alleles']].to_csv(output_file, sep='\t', index=False)


if __name__ == '__main__':
    main()

