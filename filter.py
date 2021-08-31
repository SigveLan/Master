from Bio import SeqIO
import pandas as pd
import time
from operator import itemgetter

def read_exons_to_df(file):

    # Reading to dictionary then convert to dataframe is very fast, which is why it is done here
    dictionary = {}
    ind = 0

    for seq_record in SeqIO.parse(open(file, mode='r'), 'fasta'):

        seq_record_attributes = seq_record.description.split("|")
        seq_record_attributes[5] = seq_record_attributes[5].split(':')

        # These are always present as numbers.
        for i in [3, 4, 7, 8, 13]:
            seq_record_attributes[i] = int(seq_record_attributes[i])

        # The phases can be empty, which means the same as phase zero.
        for i in [9, 10]:
            if seq_record_attributes[i] != '':
                seq_record_attributes[i] = int(seq_record_attributes[i])
            else:
                seq_record_attributes[i] = 0

        # There can be more than one coding start/stop in a given exon. This picks the closest to the edge of the exon.
        for i in [11, 12]:
            if seq_record_attributes[i] != '':
                pos = list(map(int, seq_record_attributes[i].split(';')))
                if len(pos) > 1:
                    print(seq_record_attributes)
                    differences = [abs(j - seq_record_attributes[i-4]) for j in pos]
                    seq_record_attributes[i] = pos[min(enumerate(differences), key=itemgetter(1))[0]]

                else:
                    seq_record_attributes[i] = int(seq_record_attributes[i])

        # If the sequence should be included in the dataframe - not necessary here.
        # ..Remember to add a new column in "create_df".
        #seq_record_attributes.append(seq_record.seq)

        dictionary[ind] = seq_record_attributes
        ind += 1

    return create_df(dictionary)

def create_df(dictionary):

    df = pd.DataFrame.from_dict(dictionary, orient='index', columns=["gene_id", "gene_name", "chrom", "gene_start",
                                                                     "gene_end", "transcript_ids", "exon_id",
                                                                     "exon_start", "exon_end", "phase_start",
                                                                     "phase_end", "coding_start", "coding_end",
                                                                     "strand"])
    return df


def read_SNPs_to_df(file):
    return pd.read_table(file, index_col=0)


def check_if_SNP_in_gene(SNP_pos):
    return exons[(exons['gene_start'] <= SNP_pos) & (exons['gene_end'] >= SNP_pos)]


def check_if_SNP_in_exons(gene, SNP_pos):
    return gene[(gene['exon_start'] <= SNP_pos) & (gene['exon_end'] >= SNP_pos)]


def check_if_SNP_in_coding_region(affected_exon, SNP_pos):

    if not(isinstance(affected_exon['coding_start'], int) & isinstance(affected_exon['coding_end'], int)):
        return False
    elif (affected_exon['coding_start'] <= SNP_pos) & (affected_exon['coding_end'] >= SNP_pos):
        return True
    else:
        return False


def check_SNP_location(SNP_pos):
    genes = check_if_SNP_in_gene(SNP_pos)

    if genes.empty:
        return []

    # May be multiple genes
    gene_ids = genes['gene_id'].unique()
    ids = []

    for gene_id in gene_ids:

        gene = genes.loc[genes['gene_id'] == gene_id]
        gene_name_id = [gene.iloc[0, 1], gene.iloc[0, 0]]

        affected_exons = check_if_SNP_in_exons(gene, SNP_pos)

        if len(affected_exons) == 0:
            ids.append(gene_name_id + [[], '', 'non_coding_region'])
        else:
            for i, row in affected_exons.iterrows():
                SNP_data = gene_name_id + [row['transcript_ids'], row['exon_id']]

                if check_if_SNP_in_coding_region(row, SNP_pos):
                    SNP_data.append('coding_region')

                else:
                    SNP_data.append('transcript_non_coding_region')

                ids.append(SNP_data)
    return ids


model_file = 'C:/Users/Sigve/Genome_Data/exon_model_data/exons_chrom_1.fa'
SNPs_file = 'C:/Users/Sigve/Genome_Data/SNP_data/biomart_chrom_1_test_1.tsv'

start_time = time.time()
exons = read_exons_to_df(model_file)
end_time = time.time()
print('Reading execution time = %.6f seconds' % (end_time-start_time))

SNPs = read_SNPs_to_df(SNPs_file).drop_duplicates()

start_time2 = time.time()
num_SNPs = 0

result_file_names = ['SNPs_non_coding', 'SNPs_coding']
SNP_results = [[], []]

for index, SNP in SNPs.iterrows():

    if len(SNP['Variant alleles']) != 3:
        #Skips non-SNP mutations
        continue

    location = check_SNP_location(SNP['Chromosome/scaffold position start (bp)'])

    for sublist in location:
        sublist = [index, SNP['Chromosome/scaffold name'], SNP['Chromosome/scaffold position start (bp)'], SNP['Variant alleles']] + sublist

        if sublist[8] == 'coding_region':
            SNP_results[1].append(sublist)
        else:
            SNP_results[0].append(sublist)

    num_SNPs += 1

# Write to file.
for name, result in zip(result_file_names, SNP_results):
    df = pd.DataFrame(result, columns=['variant_name', 'chrom', 'chrom_pos', 'variant_alleles', 'gene_name',
                                       'gene_id', 'transcript_ids',
                                       'exon_id', 'location'])

    df.to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/results/%s.tsv' % name, sep='\t')

end_time2 = time.time()
print('Mapping execution time = %.6f seconds' % (end_time2-start_time2))
print('Number of viable SNPs: ' + str(num_SNPs))


