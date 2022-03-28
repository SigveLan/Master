import pandas as pd

"""Writes and prints interesting results"""


def main():
    tissue_list = ['brain']#['skin', 'spleen', 'adipose_tissue', 'adrenal_gland', 'uterus', 'pancreas', 'thyroid', 'blood', 'brain', 'heart', 'kidney', 'liver', 'muscle', 'nerve']

    for tissue in tissue_list:

        input_file = 'C:/Users/Sigve/Genome_Data/results/ind_results/extensive/ind_all_chrom_sel_{0}_TEST.tsv'.format(tissue)
        output_file = 'C:/Users/Sigve/Genome_Data/results/ind_results/filtered/ind_all_chrom_sel_{0}_f_TEST.tsv'.format(tissue)

        res = pd.read_table(input_file, index_col=0)
        res['tasks_results'] = res['tasks_results'].apply(lambda x: [int(i) for i in x[1:-1].split(', ')] if x[2:-2] != 'ALL PASS' else ['ALL PASS'])
        res['solution'] = res['solution'].apply(lambda x: round(x, 3))

        ref_sol = res.at[0, 'solution']
        res['chk'] = res['tasks_results'].apply(lambda x: not all(x))

        non_nom = res[((res['solution'] != ref_sol) | res['chk']) | (res['sample_ids'] == 'REF')]
        non_nom[['sample_ids', 'solution', 'tasks_results']].to_csv(path_or_buf=output_file, sep='\t')
        print(tissue + ': ' + str(res.shape[0]))
        print(tissue + ': ' + str(non_nom.shape[0]) + ' results with damage.')
        print('-----------------------------------')


if __name__ == '__main__':
    main()
