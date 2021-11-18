import pandas as pd
from tqdm import tqdm


def SNP_filter(exons: pd.DataFrame, SNPs: pd.DataFrame) -> pd.DataFrame:
    """This script reads in SNPs and exon/gene data, then checks if the SNPs are located in an exon/gene.
       Outputs two files: one with SNPs inside exons and one with SNPs inside gene-regions
       and transcripts but not in exons."""

    def check_if_SNP_in_coding_region(affected_exon: pd.DataFrame, SNP_pos: int) -> bool:
        """Takes in an affected exon and checks if the SNP is within the coding region."""

        if not(isinstance(affected_exon['coding_start'], int) & isinstance(affected_exon['coding_end'], int)):
            return False
        elif (affected_exon['coding_start'] <= SNP_pos) & (affected_exon['coding_end'] >= SNP_pos):
            return True
        else:
            return False

    def check_SNP_location(SNP_pos: int, chrom: str) -> list:
        """Checks if a SNP is within a gene area."""

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
                        SNP_data.append('coding_sequence')

                    else:
                        SNP_data.append('transcript_non_coding')

                    ids.append(SNP_data)
        return ids

    result_list = []

    for index, SNP in tqdm(SNPs.iterrows(), total=SNPs.shape[0]):

        location = check_SNP_location(SNP['Chromosome/scaffold position start (bp)'], SNP['Chromosome/scaffold name'])

        for sublist in location:
            sublist = [index, SNP['Chromosome/scaffold name'], SNP['Chromosome/scaffold position start (bp)'], SNP['Variant alleles']] + sublist

            result_list.append(sublist)

    # Returns a dataframe of all results
    return pd.DataFrame(result_list, columns=['variant_name', 'chrom', 'chrom_pos', 'variant_alleles', 'gene_name',
                                              'gene_id', 'transcript_ids', 'exon_id', 'location'])



