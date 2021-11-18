from src.SNP_filter import SNP_filter
from src.SNP_effect_eval import SNP_effect_eval
from src.assorted_functions import SNP_sort, read_exons_to_df, add_cbm_id, df_to_tsv, split_filter_results_by_location
from src.mp_functions import parallelize_dataframe, split_filter
import time
from functools import partial


def main():

    path = 'C:/Users/Sigve/Genome_Data/'
    genome_data_file = path + 'exon_model_data/exons_chrom_all_filtered.fa'
    SNPs_file = path + 'SNP_data/SNPs_all_chrom.tsv'
    cbm_model_data_file = path + 'exon_model_data/recon-store-genes.tsv'
    output_file_names = [path + 'results/SNPs_coding.tsv',
                         path + 'results/SNPs_non_coding.tsv',
                         path + 'results/SNPs_transcript_non_coding.tsv',
                         path + 'results/SNPs_effect.tsv',
                         path + 'SNP_data/SNPs_sorted.tsv']

    write_results_to_file = True
    n_cores = 8

    start_time = time.time()

    SNPs = SNP_sort(SNPs_file)
    if write_results_to_file:
        df_to_tsv(SNPs, output_file_names[4])

    exons = read_exons_to_df(genome_data_file)
    end_time1 = time.time()
    print('Exons and SNP data loading time: %.6f seconds.' % (end_time1-start_time))

    print("Total of " + str(SNPs.shape[0]) + " SNPs divided over " + str(n_cores) + " CPU threads for mapping.")
    SNPs_in_coding = parallelize_dataframe(SNPs, partial(split_filter, partial(SNP_filter, exons)), n_cores)

    # Makes it a list of three dataframes
    SNPs_in_coding = split_filter_results_by_location(SNPs_in_coding)

    if write_results_to_file:
        for i in range(3):
            df_to_tsv(SNPs_in_coding[i], output_file_names[i])

    # Returns it to a single dataframe, now with only SNPs in coding region
    SNPs_in_coding = SNPs_in_coding[0]

    end_time2 = time.time()
    print('Mapping time: %.6f seconds\n\nNumber of SNPs in coding region: ' % (end_time2-end_time1) +
          str(SNPs_in_coding.shape[0]) + "\nChecking for changes in amino acid sequence: ")

    SNPs_effect = SNP_effect_eval(exons, SNPs_in_coding)
    end_time3 = time.time()
    print('SNP effect evaluation time: %.6f seconds' % (end_time3 - end_time2) +
          '\nNumber of affected transcripts: ' + str(SNPs_effect.shape[0]))

    SNPs_effect = add_cbm_id(cbm_model_data_file, SNPs_effect)
    if write_results_to_file:
        df_to_tsv(SNPs_effect, output_file_names[3])

    end_time4 = time.time()
    print('Total time: %.6f seconds' % (end_time4 - start_time))


if __name__ == '__main__':
    main()
