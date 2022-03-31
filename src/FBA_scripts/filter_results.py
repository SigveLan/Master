import pandas as pd

"""Writes and prints all non nominal results"""


def main():
    tissue_list = ['lung'] #['pancreas', 'spleen', 'adipose_tissue', 'adrenal_gland', 'pituitary', 'thyroid', 'blood', 'brain', 'heart', 'kidney', 'liver', 'muscle', 'nerve', 'lung']

    for tissue in tissue_list:

        input_file = 'C:/Users/Sigve/Genome_Data/results/ind_results/extensive/ind_{0}.tsv'.format(tissue)
        output_file = 'C:/Users/Sigve/Genome_Data/results/ind_results/filtered/ind_all_homozygote_{0}_f.tsv'.format(tissue)

        res = pd.read_table(input_file, index_col=0)
        res['tasks_results'] = res['tasks_results'].apply(lambda x: [int(i) for i in x[1:-1].split(', ')] if x[2:-2] != 'ALL PASS' else ['ALL PASS'])
        res['solution'] = res['solution'].apply(lambda x: round(x, 3))

        ref_sol = res.at[0, 'solution']
        res['chk'] = res['tasks_results'].apply(lambda x: not all(x))

        non_nom = res[((res['solution'] != ref_sol) | res['chk']) | (res['sample_ids'] == 'REF')]

        print(tissue + ': ' + str(res.shape[0]))
        print(tissue + ': ' + str(non_nom.shape[0] - 1) + ' results with damage.')
        print('-----------------------------------')
        if non_nom.shape[0] == 0:
            continue
        # Add gene_ids below to include them.
        non_nom[['sample_ids', 'solution', 'tasks_results']].to_csv(path_or_buf=output_file, sep='\t')


if __name__ == '__main__':
    main()
