import pandas as pd

"""Script to get SNP ids from sample ids and gene ids. Input needs sample id and gene id columns."""


def main():
    file = 'C:/Users/Sigve/Genome_Data/results/ind_results/filtered/start_stop_het_inc_ess/selected_results/{0}.tsv'
    ind_data = pd.read_table('C:/Users/Sigve/Genome_Data/results/ind_combinations/all_missense_ind_data.tsv', index_col=0)
    ind_data['gene_id'] = ind_data['gene_id'].apply(lambda x: x.split('.')[0])

    tissue_list = ['brain'] #['pancreas', 'spleen', 'adipose_tissue', 'adrenal_gland', 'pituitary', 'thyroid', 'blood', 'brain', 'heart', 'kidney', 'liver', 'muscle', 'nerve', 'lung']

    tissue_list.sort()

    for tissue in tissue_list:

        # Can be done in different ways. Either for a specific gene or whole results.
        # Can also do multiple tissues at once if input is in different tissue files

        gene = 'ENSG00000137992'
        fba_filter_res = pd.read_table(file.format(gene), index_col=0)

        #fba_filter_res = pd.read_table('C:/Users/Sigve/Genome_Data/results/ind_results/filtered/all_homozygote/subcombs/HG02078_brain_res.txt', index_col=0)


        fba_filter_res['sample_ids'] = fba_filter_res['sample_ids'].apply(lambda x: x.split(';'))
        fba_filter_res['gene_ids'].iat[0] = ''

        fba_filter_res['gene_ids'] = fba_filter_res['gene_ids'].apply(lambda x: x.split(';'))
        fba_filter_res = fba_filter_res.explode('sample_ids', ignore_index=True)

        # For selected gene, remove if to process as it in file.
        fba_filter_res['gene_ids'] = fba_filter_res['gene_ids'].apply(lambda x: [gene])

        res_list = [{'NaN': ';'}]
        for i, data in fba_filter_res.iterrows():

            if i == 0:
                continue

            res_dict = {}
            for gene in data[3]:
                gene_data = ind_data[ind_data['gene_id'] == gene]
                gene_data = ';'.join(set((gene_data[gene_data[data[0]] > 0]['variant_name'].tolist())))
                res_dict[gene] = gene_data

            res_list.append(res_dict)

        fba_filter_res['variants'] = pd.Series(res_list)
        fba_filter_res[['sample_ids', 'solution', 'tasks_results', 'variants']].to_csv(path_or_buf=file.format(gene), sep='\t')


if __name__ == '__main__':
    main()
