

import pandas as pd


"""
A simple script for producing percentages of essential genes in the different tissue models.
"""
path = 'C:/Users/Sigve/Genome_Data/results/'

tissue_list = ['spleen', 'adipose_tissue', 'adrenal_gland', 'pituitary', 'thyroid', 'blood', 'brain', 'heart',
               'kidney', 'liver', 'muscle', 'nerve', 'lung', 'skin', 'cervix_uteri', 'breast', 'blood_vessel',
               'colon', 'fallopian_tube', 'ovary', 'pancreas', 'prostate', 'salivary_gland',
               'small_intestine', 'stomach', 'testis', 'uterus', 'vagina', 'esophagus', 'bladder']

tissue_list.sort()
tissue_list = ['human1'] + tissue_list
res_list = []

for tissue in tissue_list:

    df_ess = pd.read_table(path+'model_tests/essential_genes/{0}_essential.tsv'.format(tissue), index_col=0)
    df_non_ess = pd.read_table(path + 'model_tests/essential_genes/{0}_non_essential.tsv'.format(tissue), index_col=0)

    num_ess = df_ess.shape[0]
    num_non_ess = df_non_ess.shape[0]

    res_list.append([tissue + '_'*(14-len(tissue)), num_ess, num_non_ess, num_ess+num_non_ess, round(100*(num_ess/(num_ess + num_non_ess)), 2)])


pd.DataFrame(res_list, columns=['Tissue', 'Num Essential', 'Num Nonessential', 'Total Genes', 'Percent Essential']).to_csv(path + 'model_tests/summary_essential_genes.tsv', sep='\t')



