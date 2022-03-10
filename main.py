import pandas as pd

from src.SNP_filter import SNP_filter
from src.SNP_effect_eval import SNP_effect_eval
from src.assorted_functions import SNP_sort, read_exons_to_df, add_cbm_id, df_to_tsv, split_df_results_by_column
from src.mp_functions import parallelize_dataframe, split_filter
import time
from functools import partial


def main():

    # File paths/names
    path = 'C:/Users/Sigve/Genome_Data/'
    genome_data_file = path + 'exon_model_data/exons_pc_ensembl_canonical_filtered.fa'
    SNPs_file = path + 'SNP_data/1000_genomes/result_chrom_22_snps.tsv'
    output_file_names = [path + 'results/SNPs_coding.tsv',
                         path + 'results/SNPs_non_coding.tsv',
                         path + 'results/SNPs_transcript_non_coding.tsv',
                         path + 'results/SNPs_missense.tsv',
                         path + 'results/SNPs_synonymous.tsv']

    # Settings
    write_results_to_file = True
    n_cores = 4

    start_time = time.time()

    SNPs_df = pd.read_table(SNPs_file)
    exons = read_exons_to_df(genome_data_file)

    end_time1 = time.time()
    print('Exons and SNP data loading time: %.6f seconds.' % (end_time1-start_time))

    print("Total of " + str(SNPs_df.shape[0]) + " SNPs divided over " + str(n_cores) + " CPU threads for mapping.")
    SNPs_in_coding = parallelize_dataframe(SNPs_df, partial(split_filter, partial(SNP_filter, exons)), n_cores)

    # Splits the data frame into a list of three dataframes based on location.
    # Needs to be done here due to parallelization.
    SNPs_in_coding = split_df_results_by_column(SNPs_in_coding, 'location', ['coding_sequence', 'non_coding_region', 'transcript_non_coding'])

    if write_results_to_file:
        for i in range(3):
            df_to_tsv(SNPs_in_coding[i], output_file_names[i])

    # Returns the variable to a single dataframe, now with only SNPs in coding region
    SNPs_in_coding = SNPs_in_coding[0]

    end_time2 = time.time()
    print('Mapping time: %.6f seconds\n\nNumber of SNPs in coding region: ' % (end_time2-end_time1) +
          str(SNPs_in_coding.shape[0]) + "\nChecking for changes in amino acid sequence: ")

    SNPs_effect = SNP_effect_eval(exons, SNPs_in_coding)
    end_time3 = time.time()

    # Split into synonymous and non-synonymous SNPs
    SNPs_effect = split_df_results_by_column(SNPs_effect, 'SNP_type', ['non_synonymous', 'synonymous'])
    print('SNP effect evaluation time: %.6f seconds' % (end_time3 - end_time2) +
          '\nNumber of SNPs causing amino acid change: ' + str(SNPs_effect[0].shape[0]) +
          '\nNumber of synonymous SNPs: ' + str(SNPs_effect[1].shape[0]))

    if write_results_to_file:
        df_to_tsv(SNPs_effect[0].drop(columns=['SNP_type']), output_file_names[3])
        df_to_tsv(SNPs_effect[1].drop(columns=['SNP_type', 'amino_acid_change', 'amino_acid_pos', 'score']), output_file_names[4])

    end_time4 = time.time()
    print('Total time: %.6f seconds' % (end_time4 - start_time))


if __name__ == '__main__':
    main()
