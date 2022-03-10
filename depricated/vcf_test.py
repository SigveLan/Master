import vcf
import vcf.filters
import pandas as pd
from src.assorted_functions import read_exons_to_df

"""Script for reading in VCF-file and filtering by genetic variation in genes. 
Genes have also been filtered by the model beforehand."""


vcf_reader = vcf.Reader(open('C:/Users/Sigve/Genome_Data/SNP_data/1000_genomes/ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf', 'r'))

exons = read_exons_to_df('C:/Users/Sigve/Genome_Data/exon_model_data/exons_pc_ensembl_canonical_filtered.fa')
exons = exons[exons['chrom'] == '11']

result_dict = {}
n = 0
ind = 0
for record in vcf_reader:

    if exons[(exons['gene_start'] <= record.pos) & (exons['gene_end'] >= record.pos)].empty:
        continue

    if len(record.REF) == 1 and len(record.ALT[0]) == 1:

        data = [record.POS, record.REF + '/' + str(record.ALT[0])]+[sample.gt_type for sample in record.samples]

        result_dict[ind] = data
        ind += 1


pd.DataFrame.from_dict(result_dict, orient='index', columns=['pos', 'var']+vcf_reader.samples).to_csv(
    path_or_buf='C:/Users/Sigve/Genome_Data/SNP_data/1000_genomes/result_chrom_XX.tsv', sep='\t')

