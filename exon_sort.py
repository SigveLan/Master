from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from operator import itemgetter


exon_file = 'C:/Users/sigve/Documents/Genome_Data/exon_model_data/exons_chrom_1.fa'
output_file = 'C:/Users/sigve/Documents/Genome_Data/exon_model_data/exons_chrom_1_sorted.fa'

chromosomes = [str(num) for num in range(1, 23)] + ['X', 'Y']
exon_data_by_chromosome = [[] for num in range(24)]

count_in = 0
count_out = 0

for seq_record in SeqIO.parse(open(exon_file, mode='r'), 'fasta'):

    seq_record_attributes = seq_record.description.split('|')
    exon_data_by_chromosome[chromosomes.index(seq_record_attributes[2])].append((int(seq_record_attributes[7]), int(seq_record_attributes[8]), seq_record))
    count_in += 1

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