from Bio import SeqIO
import pandas as pd

model_file = 'C:/Users/sigve/Documents/Genome_Data/biomart_queries/exons_chrom_22.fa'
#genome_file = 'C:/Users/sigve/Documents/Genome_Data/chromosome_1/Homo_sapiens.GRCh38.dna.chromosome.1.fa'
#cDNA_file = 'C:/Users/sigve/Documents/Genome_Data/cDNA_sequences/Homo_sapiens.GRCh38.cdna.all.fa'
#cds_file = 'C:/Users/sigve/Documents/Genome_Data/coding_sequences/Homo_sapiens.GRCh38.cds.all.fa'
SNPs_file = 'C:/Users/sigve/Documents/Genome_Data/SNP_data/biomarttest2.tsv'


def read_exons_to_df(file):

    df = pd.DataFrame(columns=["gene_id", "gene_name", "chrom", "gene_start", "gene_end", "transcript_ids", "exon_id",
                               "exon_start", "exon_end", "phase_start", "phase_end", "coding_start", "coding_end",
                               "strand", "sequence"])

    for seq_record in SeqIO.parse(open(file, mode='r'), 'fasta'):

        seq_record_attributes = seq_record.description.split("|")
        seq_record_attributes[5] = seq_record_attributes[5].split(':')

        for i in [3, 4, 7, 8, 9, 10, 11, 12, 13]:
            if seq_record_attributes[i] != '':
                seq_record_attributes[i] = int(seq_record_attributes[i])

        seq_record_attributes.append(seq_record.seq)
        df.loc[len(df)] = seq_record_attributes

    return df


def read_SNPs_to_df(file):
    return pd.read_table(file, index_col=0)


def check_if_SNP_in_exon(SNP_pos):
    return exons[(exons['exon_start'] <= SNP_pos) & (exons['exon_end'] >= SNP_pos)]


def get_transcripts(affected_exons):
    global exons
    t_ids = []
    print('Affected exons: ' + str(len(affected_exons)))
    for i, row in affected_exons.iterrows():
        for t_id in row['transcript_ids'][0].split(';'):
            t_ids.append(t_id)


    return t_ids


exons = read_exons_to_df(model_file)
SNPs = read_SNPs_to_df(SNPs_file)
# Transcripts is filler value.
transcripts = []
num_SNPs = 0

for index, SNP in SNPs.iterrows():

    if len(SNP['Variant alleles']) != 3:
        SNPs.drop(index)
        continue

    affected_exons = check_if_SNP_in_exon(SNP['Chromosome/scaffold position start (bp)'])

    if len(affected_exons != 0):
        transcripts = get_transcripts(affected_exons)
        for t_id in transcripts:
            print(t_id)
            print('-------------------------')
    #Check if coding change.



    num_SNPs += 1


print('Number of viable SNPs: ' + str(num_SNPs))