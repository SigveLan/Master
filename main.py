from src.SNP_filter import SNP_filter
from src.SNP_effect_eval import SNP_effect_eval
from src.assorted_functions import SNP_sort, read_exons_to_df, add_recon3d_id
import time


def main():

    start_time = time.time()

    genome_data_file = 'C:/Users/Sigve/Genome_Data/exon_model_data/exons_chrom_all_filtered.fa'
    SNPs_file = 'C:/Users/Sigve/Genome_Data/SNP_data/biomart_chrom_all_test_orig.tsv'
    recon_model_data_file = 'C:/Users/Sigve/Genome_Data/exon_model_data/recon-store-genes.tsv'

    start_time2 = time.time()
    exons = read_exons_to_df(genome_data_file)
    end_time2 = time.time()
    print('Exon loading execution time: %.6f seconds' % (end_time2-start_time2))

    SNPs = SNP_sort(SNPs_file)
    SNPs_in_coding = SNP_filter(exons, SNPs)
    SNPs_effect = SNP_effect_eval(exons, SNPs_in_coding)
    SNPs_effect = add_recon3d_id(recon_model_data_file, SNPs_effect)

    end_time = time.time()

    print('Total time: %.6f seconds' % (end_time - start_time))


if __name__ == '__main__':
    main()