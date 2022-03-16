from pysam import TabixFile, asTuple
import pandas as pd

import time


start_time = time.time()
chroms = [str(i) for i in range(1, 23)] + ['X']

missense_snps = pd.read_table('/mnt/c/Users/Sigve/Genome_Data/results/SNPs_missense.tsv', index_col=0,
                              dtype={'variation_name': str, 'chrom': str, 'chrom_pos': int, 'variant_alleles': str,
                                     'gene_name': str, 'gene_id': str, 'transcript_id': str, 'exon_id': str,
                                     'amino_acid_change': str, 'amino_acid_pos': int, 'protein_length': int,
                                     'score': float})

ind_data = []

for chrom in chroms:

    n = 0
    result_dict = {}
    _snps = missense_snps[missense_snps['chrom'] == chrom]
    vcf_in = TabixFile(
        "/mnt/c/Users/Sigve/Genome_Data/SNP_data/1000_genomes/20190312/ALL.chr{0}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz".format(chrom))

    _snps['ind_data'] = _snps['chrom_pos'].apply(lambda x: [2 if gt == '1|1' else 1 if (gt == '1|0' or gt == '0|1') else 0 for gt in next(vcf_in.fetch(chrom, x-1, x, parser=asTuple()))[9:]])

    ind_data.append(_snps)


ind_data = pd.concat(ind_data)
ind_data[list(vcf_in.header)[-1].split('\t')[9:]] = pd.DataFrame(ind_data.ind_data.to_list(), index=ind_data.index)
ind_data.drop('ind_data', axis=1, inplace=True)
ind_data.to_csv(path_or_buf='/mnt/c/Users/Sigve/Genome_Data/results/ind_combinations/test_all_chr.tsv', sep='\t')

end_time1 = time.time()
print('Variation allocation time: %.6f seconds.' % (end_time1-start_time))

