import pandas as pd

"""Script to get SNP ids from sample ids and gene ids. """


def main():

    ind_data = pd.read_table('C:/Users/Sigve/Genome_Data/results/ind_combinations/all_missense_ind_data.tsv', index_col=0)
    ind_data['gene_id'] = ind_data['gene_id'].apply(lambda x: x.split('.')[0])

    tissue_list = ['lung']#, 'spleen', 'adipose_tissue', 'adrenal_gland', 'pituitary', 'thyroid', 'blood', 'brain', 'heart', 'kidney', 'liver', 'muscle', 'nerve', 'lung']

    for tissue in tissue_list:

        fba_filter_res = pd.read_table('C:/Users/Sigve/Genome_Data/results/ind_results/filtered/ind_start_stop_{0}_f.tsv'.format(tissue), index_col=0)

        fba_filter_res['sample_ids'] = fba_filter_res['sample_ids'].apply(lambda x: x.split(';'))
        fba_filter_res['gene_ids'].iat[0] = ''

        fba_filter_res['gene_ids'] = fba_filter_res['gene_ids'].apply(lambda x: x.split(';'))
        fba_filter_res = fba_filter_res.explode('sample_ids', ignore_index=True)

        res_list = [{'NaN': ';'}]
        for i, data in fba_filter_res.iterrows():

            if i == 0:
                continue

            res_dict = {}
            for gene in data[3]:
                gene_data = ind_data[ind_data['gene_id'] == gene]
                gene_data = ';'.join(gene_data[gene_data[data[0]] == 2]['variant_name'].tolist())
                res_dict[gene] = gene_data

            res_list.append(res_dict)

        fba_filter_res['variants'] = pd.Series(res_list)
        fba_filter_res[['sample_ids', 'solution', 'tasks_results', 'variants']].to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/results/ind_results/filtered/ind_{0}_f_TEST.tsv'.format(tissue), sep='\t')


if __name__ == '__main__':
    main()
