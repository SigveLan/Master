import pandas as pd
import numpy as np

"""Simple script to write a file containing average number of genes that are essential, non essential, 
or not in tissue """

tissue_list = ['human1']#['spleen', 'adipose_tissue', 'adrenal_gland', 'pituitary', 'thyroid', 'blood', 'brain', 'heart',
               #'kidney', 'liver', 'muscle', 'nerve', 'lung', 'skin', 'cervix_uteri', 'breast', 'blood_vessel',
               #'colon', 'fallopian_tube', 'ovary', 'pancreas', 'prostate', 'salivary_gland',
               #'small_intestine', 'stomach', 'testis', 'uterus', 'vagina', 'esophagus', 'bladder']

tissue_list.sort()


path = 'C:/Users/Sigve/Genome_Data/results/'
input_file = path + 'ind_combinations/start_stop_comb_het.tsv'


ind_data = pd.read_table(input_file, index_col=0)
ind_data['gene_ids'] = ind_data['gene_ids'].apply(lambda x: x[2:-2].split(';'))

results = []

for tissue in tissue_list:

    genes_non_ess = \
    pd.read_table(path + 'model_tests/essential_genes/{0}_non_essential.tsv'.format(tissue), index_col=0)[
        'gene_ids'].tolist()
    genes_ess = \
    pd.read_table(path + 'model_tests/essential_genes/{0}_essential.tsv'.format(tissue), index_col=0)[
        'gene_ids'].tolist()

    ind_data[['sample_ids', 'gene_ids']].copy()

    ind_data['total'] = ind_data['gene_ids'].apply(len)
    ind_data['ess'] = ind_data['gene_ids'].apply(lambda x: len(list(set(x).intersection(genes_ess))))
    ind_data['non_ess'] = ind_data['gene_ids'].apply(lambda x: len(list(set(x).intersection(genes_non_ess))))
    ind_data['not_in_tissue'] = ind_data['total'] - ind_data['ess'] - ind_data['non_ess']

    res = [tissue, np.average(ind_data['total'].tolist()),
                            np.average(ind_data['non_ess'].tolist()),
                            np.average(ind_data['ess'].tolist()),
                            np.average(ind_data['not_in_tissue'].tolist())]

    results.append(res + [(res[2]/res[1])*100, (res[3]/res[1])*100, (res[4]/res[1])*100])

results = pd.DataFrame(results, columns=['tissue', 'total', 'non_ess', 'ess', 'not_in_tissue', 'percent_non_ess', 'percent_ess', 'percent_not_in_tissue'])
results.to_csv(path_or_buf=path+'ind_results/start_stop/human_1_all_start_stop_het_genes_sample.tsv', sep='\t')