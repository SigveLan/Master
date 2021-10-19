import pandas as pd
import time
from src.assorted_functions import read_exons_to_df


def check_if_SNP_in_gene(SNP_pos: int, chrom: str, exons: pd.DataFrame) -> str:
    """Checks if a SNP is within a gene region, and if it is, returns that gene"""

    gene_df = exons[(exons['gene_start'] <= SNP_pos) & (exons['gene_end'] >= SNP_pos)]
    gene_df = gene_df[gene_df['chrom'] == chrom]

    if gene_df.empty:
        return 'Not in gene'
    else:
        return gene_df['gene_name'].iloc[0]


def main ():

    model_file = 'C:/Users/Sigve/Genome_Data/exon_model_data/exons_chrom_all_filtered.fa'
    SNPs_file = 'C:/Users/Sigve/Genome_Data/SNP_data/biomart_chrom_all_test.tsv'

    # Timer for reading in exon data file, not necessary for function
    start_time = time.time()

    exons = read_exons_to_df(model_file)
    SNP_df = pd.read_table(SNPs_file)

    end_time = time.time()
    print('Reading execution time = %.6f seconds' % (end_time-start_time))

    start_time2 = time.time()

    num_SNPs = 0

    result_files_names = ['SNPs_non_coding', 'SNPs_coding']
    SNP_results = [[], []]

    SNP_df['gene_name'] = SNP_df.apply(check_if_SNP_in_gene(SNP_df['Chromosome/scaffold position start (bp)'], SNP_df['Chromosome/scaffold name'], exons))
    SNP_df = SNP_df[SNP_df['gene_name'] != 'Not in gene']


    print(SNP_df.head())
    exit()




    for name, result in zip(result_files_names, SNP_results):
        df = pd.DataFrame(result, columns=['variant_name', 'chrom', 'chrom_pos', 'variant_alleles', 'gene_name',
                                           'gene_id', 'transcript_ids',
                                           'exon_id', 'location'])
        # Change the "%" to new format
        df.to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/results/%s.tsv' % name, sep='\t')

    end_time2 = time.time()
    print('Mapping execution time = %.6f seconds' % (end_time2 - start_time2))
    print('Number of viable SNPs: ' + str(num_SNPs))


if __name__ == '__main__':
    main()

# >ENSG00000176022.7|B3GALT6|1|1232237|1235041|ENST00000379198.5|ENSE00001480062|1232237|1235041|-1|-1|1232279|1233268|1
"""columns=["gene_id", "gene_name", "chrom", "gene_start",
                                                                     "gene_end", "transcript_ids", "exon_id",
                                                                     "exon_start", "exon_end", "phase_start",
                                                                     "phase_end", "coding_start", "coding_end",
                                                                     "strand", "sequence"])"""