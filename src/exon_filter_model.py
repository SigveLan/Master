from Bio import SeqIO
from operator import itemgetter
import pandas as pd


"""Filters for genes in the model. Also sorts the exons in the datafile from Ensembl."""

# List of exons and sequence data
exon_file = 'C:/Users/Sigve/Genome_Data/exon_model_data/exons_pc_ensembl_canonical.fa'

# Model gene list
model_file = 'C:/Users/Sigve/Genome_Data/Human1/Human1_GEM/genes.tsv'

# Output file
output_file = 'C:/Users/Sigve/Genome_Data/exon_model_data/exons_pc_ensembl_canonical_filtered.fa'


model_gene_ids = pd.read_table(model_file)['genes'].tolist()

chromosomes = [str(num) for num in range(1, 23)] + ['X', 'Y', 'MT']
exon_data_by_chromosome = [[] for num in range(25)]

count_in = 0
count_out = 0

for seq_record in SeqIO.parse(open(exon_file, mode='r'), 'fasta'):

    seq_record_attributes = seq_record.description.split('|')
    count_in += 1
    gene_id = seq_record_attributes[0].split('.')[0]

    if gene_id not in model_gene_ids:
        continue
    exon_data_by_chromosome[chromosomes.index(seq_record_attributes[2])].append((int(seq_record_attributes[7]), int(seq_record_attributes[8]), seq_record))


with open(output_file, 'w') as f_out:

    for chromosome in exon_data_by_chromosome:
        chromosome.sort(key=itemgetter(0, 1))

        for location in chromosome:
            r = SeqIO.write(location[2], f_out, 'fasta')

            if r != 1:
                print('Error while writing sequence:  ' + location[1].id)
            count_out += 1

print("Records in: " + str(count_in))
print("Records out: " + str(count_out))