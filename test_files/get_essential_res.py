import pandas as pd


"""
A simple script for comparison of essential gene FBA runs in different tissues. 
Can probably be used for other stuff as well.

"""
path = 'C:/Users/Sigve/Genome_Data/results/'
genes = pd.read_table(path + 'model_tests/essential_genes/liver_essential.tsv')['gene_ids'].tolist()

#genes = 'ENSG00000140740;ENSG00000164405;ENSG00000169021;ENSG00000127540'.split(';')

res_list = []

for gene in genes:

    tissue_list = ['heart', 'blood', 'brain', 'muscle', 'lung', 'liver']

    for tissue in tissue_list:

        df = pd.read_table(path+'model_tests/ess_results/essential_genes_FBA_res_{0}.tsv'.format(tissue), index_col=0)

        df = df[df['gene_ids'] == gene]
        if df.shape[0] == 0:
            continue
        res_list.append([gene, tissue] + [df.solution.iat[0], df.tasks_results.iat[0]])
    res_list.append(['HurDurBreakhere', 'Skip_line', 0, ['____________________________________________________']])

pd.DataFrame(res_list, columns=['gene_ids', 'tissue', 'solution', 'tasks_results']).to_csv(path+'model_tests/ess_results/extracted_results/{0}.tsv'.format(gene), sep='\t')