import time

import cobra.io
import pandas as pd
from functools import partial

from src.mp_functions import combinations_subset, parallelize_dataframe, knockout_FBA_simple

"""Script for simple FBA. Be sure to change file in-/out-put names!"""


def main():

    start_time = time.time()

    tissue_list = ['liver']#['pancreas', 'spleen', 'adipose_tissue', 'adrenal_gland', 'pituitary', 'thyroid', 'blood', 'brain', 'heart', 'kidney', 'liver', 'muscle', 'nerve', 'lung']
    #tissue_list.sort()

    ind_data = pd.read_table('C:/Users/Sigve/Genome_Data/results/ind_combinations/HG03385_liver_combinations.tsv', index_col=0)
    ind_data['gene_ids'] = ind_data['gene_ids'].apply(lambda x: x[2:-2].split(';'))

    # Just for cleanup, not necessary for functions
    ind_data['sample_ids'] = ind_data['sample_ids'].apply(lambda x: ';'.join(x[2:-2].split("', '")))

    for tissue in tissue_list:

        # Reads in non essential genes list.
        # This list includes genes not essential for tasks, otherwise cobra get essential can be used for comparison.
        genes = pd.read_table('C:/Users/Sigve/Genome_Data/results/model_tests/essential_genes/{0}_non_essential.tsv'.format(tissue), index_col=0)['gene_ids'].tolist()

        # Load the model
        model = cobra.io.read_sbml_model('C:/Users/Sigve/Genome_Data/Human1/Human1_GEM/GTEx/{0}.xml'.format(tissue))

        # Use all genes instead of non essential genes; needed because not all tissues have the same genes.
        # genes= [gene.id for gene in model_list[0].genes]

        results = ind_data.copy()
        results['gene_ids'] = results['gene_ids'].apply(lambda x: list(set(x).intersection(genes)))
        results = results[results['gene_ids'].map(lambda x: len(x)) > 0]

        # Actual FBA
        results = parallelize_dataframe(results, partial(combinations_subset, partial(knockout_FBA_simple, model)), n_cores=12)

        results.loc[-1] = ['REF', [], knockout_FBA_simple(model, [])]

        # Put REF on top
        results.index = results.index + 1  # shifting index
        results.sort_index(inplace=True)

        results['gene_ids'] = results['gene_ids'].apply(';'.join)
        results['solution'] = results['results'].apply(lambda x: round(x, 3))

        results.reset_index(inplace=True, drop=True)
        results[['sample_ids', 'gene_ids', 'solution']].to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/results/ind_results/extensive/HG03385_{0}_res.tsv'.format(tissue), sep='\t')

    end_time = time.time()

    print('Total time: %.6f seconds' % (end_time - start_time))


if __name__ == '__main__':
    main()
