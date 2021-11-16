from src.SNP_filter import SNP_filter
from src.SNP_effect_eval import SNP_effect_eval
from src.assorted_functions import SNP_sort, read_exons_to_df, add_cbm_id
import pandas as pd

genome_data_file = 'C:/Users/Sigve/Genome_Data/exon_model_data/exons_chrom_all_filtered.fa'
SNPs_file = 'C:/Users/Sigve/Genome_Data/SNP_data/SNPs_all_chrom.tsv'
cbm_model_data_file = 'C:/Users/Sigve/Genome_Data/exon_model_data/recon-store-genes.tsv'

SNPs_coding = pd.read_table('C:/Users/Sigve/Genome_Data/results/SNPs_coding.tsv', index_col=0)

SNPs_coding['transcript_ids'] = SNPs_coding['transcript_ids'].apply(lambda x: x[2:-2].split("', '"))

exons = read_exons_to_df(genome_data_file)

res = SNP_effect_eval(exons, SNPs_coding)


