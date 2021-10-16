import pandas as pd
import time
from read_exons_to_df import read_exons_to_df
"""
This script reads in SNPs and exon data, then checks if the SNPs are located in an exon.
Outputs two files: one with SNPs inside exons and one with SNPs inside genes but not in exons.
"""
def read_SNPs_to_df(file):
    """Reads in SNP data"""
    return pd.read_table(file, index_col=0)


def check_if_SNP_in_gene(SNP_pos):
    """Checks if a SNP is within a gene region, and if it is, returns that gene"""
    return exons[(exons['gene_start'] <= SNP_pos) & (exons['gene_end'] >= SNP_pos)]


def check_if_SNP_in_exons(gene, SNP_pos):
    """Takes in all exons of a gene and checks if a SNP is within one or more, returns affected exons"""
    return gene[(gene['exon_start'] <= SNP_pos) & (gene['exon_end'] >= SNP_pos)]


def check_if_SNP_in_coding_region(affected_exon, SNP_pos):
    """Takes in an affected exon and checks if the SNP is within the coding region"""
    if not(isinstance(affected_exon['coding_start'], int) & isinstance(affected_exon['coding_end'], int)):
        return False
    elif (affected_exon['coding_start'] <= SNP_pos) & (affected_exon['coding_end'] >= SNP_pos):
        return True
    else:
        return False


def check_SNP_location(SNP_pos):
    """Checks if a SNP is within an area of interest"""
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

model_file = 'C:/Users/Sigve/Genome_Data/exon_model_data/exons_chrom_all_filtered.fa'
SNPs_file = 'C:/Users/Sigve/Genome_Data/SNP_data/biomart_chrom_all_test.tsv'

# Timer for reading in exon data file, not necessary for function
start_time = time.time()

exons = read_exons_to_df(model_file)

end_time = time.time()
print('Reading execution time = %.6f seconds' % (end_time-start_time))

SNPs = read_SNPs_to_df(SNPs_file).drop_duplicates()

start_time2 = time.time()
num_SNPs = 0

result_files_names = ['SNPs_non_coding', 'SNPs_coding']
SNP_results = [[], []]

for index, SNP in SNPs.iterrows():

    location = check_SNP_location(SNP['Chromosome/scaffold position start (bp)'])

    for sublist in location:
        sublist = [index, SNP['Chromosome/scaffold name'], SNP['Chromosome/scaffold position start (bp)'], SNP['Variant alleles']] + sublist

        if sublist[8] == 'coding_region':
            SNP_results[1].append(sublist)
        else:
            SNP_results[0].append(sublist)

    num_SNPs += 1

# Write results to file.
for name, result in zip(result_files_names, SNP_results):
    df = pd.DataFrame(result, columns=['variant_name', 'chrom', 'chrom_pos', 'variant_alleles', 'gene_name',
                                       'gene_id', 'transcript_ids',
                                       'exon_id', 'location'])
    # Change the "%" to new format
    df.to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/results/%s.tsv' % name, sep='\t')

end_time2 = time.time()
print('Mapping execution time = %.6f seconds' % (end_time2-start_time2))
print('Number of viable SNPs: ' + str(num_SNPs))


