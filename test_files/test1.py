from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

# A file used to test various code

'''
genome_file = 'C:/Users/sigve/Documents/Genome_Data/chromosome_22/Homo_sapiens.GRCh38.dna.chromosome.22.fa'
cDNA_file = 'C:/Users/sigve/Documents/Genome_Data/cDNA_sequences/Homo_sapiens.GRCh38.cdna.all.fa'
cds_file = 'C:/Users/sigve/Documents/Genome_Data/coding_sequences/Homo_sapiens.GRCh38.cds.all.fa'
model_file = 'C:/Users/sigve/Documents/Genome_Data/model/model.fa'
out_file = 'C:/Users/sigve/Documents/Genome_Data/SNP_data/rsids.txt'
'''

'''
def check_pos(pos):
    for seq_record in SeqIO.parse(open(genome_file, mode='r'), 'fasta'):
        chromosome22 = seq_record.seq
        print(chromosome22[pos-2:pos+1])
        print(chromosome22[pos-1])
check_pos(10511507)
'''
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

'''
def identify_codon(pos, start_phase, end_phase, strand, seq):
    if (1 < pos) & (pos < len(seq) - 2):
        print('...............................')
        reverse = False

        if strand == '-1':
            seq = seq[::-1]
            start_phase = end_phase
            reverse = True
            pos = abs(pos - len(seq) + 1)
            print(str(pos) + '----' + str(len(seq)))

        if start_phase != '':
            offset = abs(int(start_phase) - 3*(not reverse))
            print('Offset: ' + str(offset))
            seq = seq[offset:]
            pos = pos - offset
        print(seq)
        return seq[pos - pos%3:pos - pos%3 + 3]
    return ''
'''
def identify_codon(pos, start_phase, seq, var):

    var = var.split('/')

    if (1 < pos) & (pos < len(seq) - 2):

        if start_phase != '':
            offset = abs(int(start_phase)-3)
            pos = pos - offset
            seq = seq[offset:]

        var_seq = seq[0:pos] + var[1] + seq[pos + 1:]
        return [seq[pos - pos % 3:pos - pos % 3 + 3], var_seq[pos - pos % 3:pos - pos % 3 + 3]]




print(identify_codon(6, '', 'ATGATTATCCTG', 'A/G'))
print(identify_codon(5, '1', 'TGATTATCCTG', 'A/G'))
print(identify_codon(4, '2',  'GATTATCCTG', 'A/G'))

print(identify_codon(7, '', 'ATGATTATCCTG', 'A/G'))
print(identify_codon(6, '1', 'TGATTATCCTG', 'A/G'))
print(identify_codon(5, '2',  'GATTATCCTG', 'A/G'))

print(identify_codon(8, '', 'ATGATTATCCTG', 'A/G'))
print(identify_codon(7, '1', 'TGATTATCCTG', 'A/G'))
print(identify_codon(6, '2',  'GATTATCCTG', 'A/G'))

'''
exon = get_exon_by_id(data['exon_id'])
codons = identify_codon(data['chrom_pos'] - exon.iloc[0, 7], exon.iloc[0, 9], exon.iloc[0]['sequence'],
                        data['variant_alleles'])

if not codons:
    ancestral_AA.append('')
    var_AA.append('')
    change_score.append('')
    continue

translate_result = translate_codon(codons)
ancestral_AA.append(translate_result[0])
var_AA.append(translate_result[1])
change_score.append(compare_and_score_AA(translate_result))

'''