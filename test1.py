from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd


genome_file = 'C:/Users/sigve/Documents/Genome_Data/chromosome_22/Homo_sapiens.GRCh38.dna.chromosome.22.fa'
cDNA_file = 'C:/Users/sigve/Documents/Genome_Data/cDNA_sequences/Homo_sapiens.GRCh38.cdna.all.fa'
cds_file = 'C:/Users/sigve/Documents/Genome_Data/coding_sequences/Homo_sapiens.GRCh38.cds.all.fa'
model_file = 'C:/Users/sigve/Documents/Genome_Data/model/model.fa'
out_file = 'C:/Users/sigve/Documents/Genome_Data/SNP_data/rsids.txt'



def check_pos(pos):
    for seq_record in SeqIO.parse(open(genome_file, mode='r'), 'fasta'):
        chromosome22 = seq_record.seq
        print(chromosome22[pos-2:pos+1])
        print(chromosome22[pos-1])
check_pos(10511507)

'''
SNP_data = pd.read_table('C:/Users/sigve/Documents/Genome_Data/SNP_data/EirikUKbioSNPPhe.txt')
rsids = SNP_data.query('chrom==1')['rsids']

SNP_data_2 = pd.read_table('C:/Users/sigve/Documents/Genome_Data/SNP_data/homo_sapiens_somatic.gvf', skiprows=32, nrows=1000)
print(SNP_data_2.head(20))
'''

'''
with open(out_file, 'w') as f_out:


    count = 0
    while count < 500:
        f_out.write(rsids.iloc[count] + "\n")
        count += 1
'''