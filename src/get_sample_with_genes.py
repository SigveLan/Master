import pandas as pd


def main():
    """Get samples with specific genes"""

    tissue = 'pancreas'
    target_genes = ['ENSG00000125246']

    fba_filter_res = pd.read_table('C:/Users/Sigve/Genome_Data/results/ind_results/filtered/start_stop_het/ind_results_f_{0}.tsv'.format(tissue), index_col=0)

    fba_filter_res['gene_ids'].iat[0] = ''

    fba_filter_res['gene_ids'] = fba_filter_res['gene_ids'].apply(lambda x: x.split(';'))
    fba_filter_res = fba_filter_res.explode('sample_ids', ignore_index=True)

    ref = fba_filter_res.iloc[[0]]

    fba_filter_res = fba_filter_res[fba_filter_res['gene_ids'].map(lambda x: all([gene in x for gene in target_genes]))]
    fba_filter_res['gene_ids'] = fba_filter_res['gene_ids'].apply(lambda x: ';'.join(x))

    # To include only target genes in output
    #fba_filter_res['gene_ids'] = target_genes[0]

    fba_filter_res = pd.concat([ref, fba_filter_res])
    print(fba_filter_res)

    fba_filter_res.to_csv(
        path_or_buf='C:/Users/Sigve/Genome_Data/results/ind_results/filtered/start_stop_het/ind_results_f_{0}_selected.tsv'.format(
            tissue), sep='\t')


if __name__ == '__main__':
    main()
