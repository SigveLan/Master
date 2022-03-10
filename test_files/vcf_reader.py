from pysam import TabixFile, asTuple
import pandas as pd
from src.assorted_functions import read_exons_to_df
import time

"""Script for reading in VCF-file and filtering by genetic variation in genes. 
Genes have also been filtered by the model beforehand.  Note: this script must be run on Linux"""


def write_data(data, n, chrom, header) -> None:
    pd.DataFrame.from_dict(data, orient='index',
                           columns=['pos', 'var'] + list(header)[-1].split('\t')[9:]) \
        .to_csv(path_or_buf='/mnt/c/Users/Sigve/Genome_Data/SNP_data/1000_genomes/results_20190312/chrom_{0}_filtered_{1}.tsv'.format(chrom, n), sep='\t')


start_time = time.time()
chroms = [str(i) for i in range(1, 23)] + ['X']
exons = read_exons_to_df('/mnt/c/Users/Sigve/Genome_Data/exon_model_data/exons_pc_ensembl_canonical_filtered.fa')

ind = 0
for chrom in chroms:

    n = 1
    result_dict = {}
    _exons = exons[exons['chrom'] == chrom]
    regions = zip(_exons['gene_start'].unique().tolist(), _exons['gene_end'].unique().tolist())

    vcf_in = TabixFile(
        "/mnt/c/Users/Sigve/Genome_Data/SNP_data/1000_genomes/20190312/ALL.chr{0}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz".format(chrom))

    for start, end in regions:
        for rec in vcf_in.fetch(chrom, start-1, end+1, parser=asTuple()):

            if rec[7].split(';')[-2] != 'VT=SNP':
                continue

            data = [rec[1], rec[3] + '/' + rec[4]] + \
                   [2 if gt == '1|1' else 1 if (gt == '1|0' or gt == '0|1') else 0 for gt in rec[9:]]

            result_dict[ind] = data
            ind += 1

            if len(result_dict) == 100000:
                write_data(result_dict, n, chrom, vcf_in.header)
                result_dict = {}
                n += 1

    write_data(result_dict, n, chrom, vcf_in.header)

end_time1 = time.time()
print('Variation filter time: %.6f seconds.' % (end_time1-start_time))
