import pandas as pd

"""This script reads in the results from 'vcf_reader' and combines the SNPs into single file.
Only necessary if multiple files are produced. Will work with multiple files for each chromosome."""


def main():

    chroms = [str(i) for i in range(1, 23)] + ['X']
    res_list = []
    for chrom in chroms:
        n = 1
        while n>0:
            try:
                snp_data = pd.read_table('C:/Users/Sigve/Genome_Data/SNP_data/1000_genomes/results_20190312/SNPs_chrom_{0}_{1}.tsv'.format(chrom, n), index_col=0)
                snp_data['snp_id'] = snp_data['pos'].apply(lambda x: chrom + ':' + str(x))
                snp_data['chr'] = chrom
                res_list.append(snp_data[['snp_id', 'chr', 'pos', 'var']])
                n += 1
            except FileNotFoundError:
                break

    pd.concat(res_list, ignore_index=True).set_index('snp_id').to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/SNP_data/1000_genomes/SNP_results_all_chrom.tsv', sep='\t')


if __name__ == '__main__':
    main()
