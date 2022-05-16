from pysam import TabixFile, asTuple
import pandas as pd
from src.assorted_functions import read_exons_to_df
import time

"""Script for reading in VCF files and filtering by genetic variation in genes. 
Genes should also have been filtered by the model beforehand. NB: pysam is not available for windows."""


def write_data(data, n, chrom) -> None:
    """Writes collected SNPs to file"""
    pd.DataFrame.from_dict(data, orient='index', columns=['pos', 'var']) \
        .to_csv(path_or_buf='/mnt/c/Users/Sigve/Genome_Data/SNP_data/1000_genomes/results_20190312/SNPs_chrom_{0}_{1}.tsv'.format(chrom, n), sep='\t')


start_time = time.time()
chroms = [str(i) for i in range(1, 23)] + ['X']
exons = read_exons_to_df('/mnt/c/Users/Sigve/Genome_Data/exon_model_data/exons_pc_ensembl_canonical_filtered.fa')

# Multiple files per chromosomes depending on number of entries.
# This will reduce memory requirements with large data amounts.
multiple_files = False

ind = 0
for chrom in chroms:

    n = 1
    result_dict = {}
    _exons = exons[exons['chrom'] == chrom]

    # Could add en extended range upstream and/or downstream of the gene.
    regions = zip(_exons['gene_start'].unique().tolist(), _exons['gene_end'].unique().tolist())

    vcf_in = TabixFile(
        "/mnt/c/Users/Sigve/Genome_Data/SNP_data/1000_genomes/20190312/ALL.chr{0}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz".format(chrom))

    for start, end in regions:
        for rec in vcf_in.fetch(chrom, start-1, end+1, parser=asTuple()):

            # Only SNPs are registered
            if rec[7].split(';')[-2] != 'VT=SNP':
                continue

            data = [rec[1], rec[3] + '/' + rec[4]]

            result_dict[ind] = data
            ind += 1

            # How many entries per output file
            if multiple_files and len(result_dict) == 100000:
                write_data(result_dict, n, chrom)
                result_dict = {}
                n += 1

    write_data(result_dict, n, chrom)

end_time1 = time.time()
print('Variation filter time: %.6f seconds.' % (end_time1-start_time))
