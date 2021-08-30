from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

genome_file = 'C:/Users/sigve/Documents/Genome_Data/chromosome_1/Homo_sapiens.GRCh38.dna.chromosome.1.fa'
cDNA_file = 'C:/Users/sigve/Documents/Genome_Data/cDNA_sequences/Homo_sapiens.GRCh38.cdna.all.fa'
cds_file = 'C:/Users/sigve/Documents/Genome_Data/coding_sequences/Homo_sapiens.GRCh38.cds.all.fa'
model_file = 'C:/Users/sigve/Documents/Genome_Data/model/model.fa'


chromosome1 = ''
for seq_record in SeqIO.parse(open(genome_file, mode='r'), 'fasta'):
    chromosome1 = seq_record.seq


with open(model_file, 'w') as f_out:
    for seq_record in SeqIO.parse(open(cds_file, mode='r'), 'fasta'):

        seq_record.attributes = seq_record.description.split()

        location_info = seq_record.attributes[2].split(':')
        chromosome = location_info[2]

        if chromosome != '1':
            continue

        if seq_record.attributes[4].split(':')[1] != 'protein_coding':
            continue
        print(seq_record.attributes)

        sequence_start = int(location_info[3])
        sequence_stop = int(location_info[4])

        if chromosome1[sequence_start-1:sequence_stop] == seq_record.seq:
            print("They are equal!")
        else:
            print('They are not equal :(')
            print(len(chromosome1[sequence_start-1:sequence_stop]))
            print(len(seq_record.seq))

            print(chromosome1[sequence_start-1:sequence_start+9])
            print(seq_record.seq[:10])

            print(chromosome1[sequence_stop-10:sequence_stop])
            print(seq_record.seq[-10:])

            print(chromosome1[sequence_start - 1:sequence_stop])
            print(seq_record.seq)

        #print('SequenceID = ' + seq_record.id)
        #print('Description = ' + seq_record.description + '\n')

        #remove -id from description
        seq_record.description = ' '.join(seq_record.attributes[1:])

        # write new fasta file
        r = SeqIO.write(seq_record, f_out, 'fasta')


        if r != 1:
            print('Error while writing sequence:  ' + seq_record.id)
        break


