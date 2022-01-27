from src.assorted_functions import read_exons_to_df
from collections import Counter

path = 'C:/Users/Sigve/Genome_Data/'
genome_data_file = path + 'exon_model_data/exons_chrom_all_filtered.fa'

exons = read_exons_to_df(genome_data_file)

exons['gene_id'] = exons['gene_id'].apply(lambda x: x.split('.')[0])
id_list = exons['gene_id'].tolist()

#print(Counter(id_list).keys()) # equals to list(set(words))
#print(Counter(id_list).values()) # counts the elements' frequency

print(len(set(id_list)))


