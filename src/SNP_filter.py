import pandas as pd
import time
from tqdm import tqdm


def SNP_filter(exons: pd.DataFrame, SNPs: pd.DataFrame, write_results_to_file=True) -> pd.DataFrame:
    """This script reads in SNPs and exon data, then checks if the SNPs are located in an exon.
       Outputs two files: one with SNPs inside exons and one with SNPs inside genes but not in exons."""

    def check_if_SNP_in_coding_region(affected_exon: pd.DataFrame, SNP_pos: int) -> bool:
        """Takes in an affected exon and checks if the SNP is within the coding region"""
        if not(isinstance(affected_exon['coding_start'], int) & isinstance(affected_exon['coding_end'], int)):
            return False
        elif (affected_exon['coding_start'] <= SNP_pos) & (affected_exon['coding_end'] >= SNP_pos):
            return True
        else:
            return False

    def check_SNP_location(SNP_pos: int, chrom: str) -> list:
        """Checks if a SNP is within an area of interest"""

        genes = exons[(exons['gene_start'] <= SNP_pos) & (exons['gene_end'] >= SNP_pos)]
        genes = genes.loc[genes['chrom'] == str(chrom)]

        if genes.empty:
            return []

        # May be multiple genes
        gene_ids = genes['gene_id'].unique()
        ids = []

        for gene_id in gene_ids:

            gene = genes.loc[genes['gene_id'] == gene_id]
            gene_name_id = [gene.iloc[0, 1], gene.iloc[0, 0]]

            affected_exons = gene[(gene['exon_start'] <= SNP_pos) & (gene['exon_end'] >= SNP_pos)]

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

    def result_to_df(result_list: list) -> pd.DataFrame:
        return pd.DataFrame(result_list, columns=['variant_name', 'chrom', 'chrom_pos', 'variant_alleles', 'gene_name',
                                                   'gene_id', 'transcript_ids', 'exon_id', 'location'])

    start_time2 = time.time()

    num_SNPs = 0
    SNP_results = [[], []]

    for index, SNP in tqdm(SNPs.iterrows(), total=SNPs.shape[0]):

        location = check_SNP_location(SNP['Chromosome/scaffold position start (bp)'], SNP['Chromosome/scaffold name'])

        for sublist in location:
            sublist = [index, SNP['Chromosome/scaffold name'], SNP['Chromosome/scaffold position start (bp)'], SNP['Variant alleles']] + sublist

            if sublist[8] == 'coding_region':
                SNP_results[1].append(sublist)
                num_SNPs += 1
            else:
                SNP_results[0].append(sublist)

    SNP_results_as_dfs = [result_to_df(lst) for lst in SNP_results]

    if write_results_to_file:
        for name, result in zip(['SNPs_non_coding', 'SNPs_coding'], SNP_results_as_dfs):
            result.to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/results/{!s}.tsv'.format(name), sep='\t')

    end_time2 = time.time()
    print('\nMapping execution time: %.6f seconds' % (end_time2-start_time2))
    print('Number of SNPs in coding region: ' + str(num_SNPs))

    # Only return the in coding region results.
    return SNP_results_as_dfs[1]


