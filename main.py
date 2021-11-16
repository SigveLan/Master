from src.SNP_filter import SNP_filter
from src.SNP_effect_eval import SNP_effect_eval
from src.assorted_functions import SNP_sort, read_exons_to_df, add_cbm_id, parallelize_dataframe, split_filter, df_to_tsv
import time
from functools import partial


def main():

    start_time = time.time()

    genome_data_file = 'C:/Users/Sigve/Genome_Data/exon_model_data/exons_chrom_all_filtered.fa'
    SNPs_file = 'C:/Users/Sigve/Genome_Data/SNP_data/SNPs_all_chrom.tsv'
    cbm_model_data_file = 'C:/Users/Sigve/Genome_Data/exon_model_data/recon-store-genes.tsv'
    write_results_to_file = True

    SNPs = SNP_sort(SNPs_file)
    exons = read_exons_to_df(genome_data_file)
    end_time1 = time.time()
    print('Exon loading time: %.6f seconds.' % (end_time1-start_time))

    n_cores = 12
    print("Total of " + str(SNPs.shape[0]) + " SNPs divided over " + str(n_cores) + " CPU threads for mapping.")
    SNPs_in_coding = parallelize_dataframe(SNPs, partial(split_filter, partial(SNP_filter, exons)), n_cores)

    SNPs_in_coding.reset_index(drop=True, inplace=True)
    if write_results_to_file:
        df_to_tsv(SNPs_in_coding, 'C:/Users/Sigve/Genome_Data/results/SNPs_coding.tsv')

    end_time2 = time.time()
    print('\nMapping time: %.6f seconds.' % (end_time2-end_time1))
    print('Number of SNPs in coding region: ' + str(SNPs_in_coding.shape[0]) + "\n")

    SNPs_effect = SNP_effect_eval(exons, SNPs_in_coding)
    if write_results_to_file:
        df_to_tsv(SNPs_effect, 'C:/Users/Sigve/Genome_Data/results/SNPs_effect.tsv')
    end_time3 = time.time()
    print('SNP effect evaluation time: %.6f seconds.' % (end_time3 - end_time2))

    SNPs_effect = add_cbm_id(cbm_model_data_file, SNPs_effect)
    if write_results_to_file:
        df_to_tsv(SNPs_effect, 'C:/Users/Sigve/Genome_Data/results/SNPs_effect_mod.tsv')

    end_time4 = time.time()
    print('Total time: %.6f seconds' % (end_time4 - start_time))


if __name__ == '__main__':
    main()