import time
import numpy as np
import pandas as pd
from functools import partial

from src.mp_functions import parallelize_dataframe, knockout_FBA_w_tasks
from src.FBA_scripts.met_task_functions import read_tasks, constrain_model

"""Script for FBA with tasks. Prepares each task then divides the  Be sure to change file in- and output names!"""


def main():

    start_time = time.time()

    tissue_list = ['spleen', 'adipose_tissue', 'adrenal_gland', 'pituitary', 'thyroid', 'blood', 'brain', 'heart', 'kidney', 'liver', 'muscle', 'nerve', 'lung']
    tissue_list.sort()

    ind_data = pd.read_table('C:/Users/Sigve/Genome_Data/results/ind_combinations/start_stop_comb_het.tsv', index_col=0)
    ind_data['gene_ids'] = ind_data['gene_ids'].apply(lambda x: x[2:-2].split(';'))

    # Just for cleanup, not necessary for functions
    ind_data['sample_ids'] = ind_data['sample_ids'].apply(lambda x: ';'.join(x[2:-2].split("', '")))

    for tissue in tissue_list:

        # Reads in non essential genes list.
        # This list includes genes not essential for tasks, otherwise cobra get essential can be used for comparison.
        #genes = pd.read_table('C:/Users/Sigve/Genome_Data/results/model_tests/essential_genes/{0}_non_essential.tsv'.format(tissue), index_col=0)['gene_ids'].tolist()

        # Load models. Multiple instances are needed to do both regular and task FBA.
        model_list = constrain_model('C:/Users/Sigve/Genome_Data/Human1/Human1_GEM/GTEx/{0}.xml'.format(tissue), ALLMETSIN_OUT=True)

        # Use all genes instead; needed because not all tissues have the same genes.
        genes = [gene.id for gene in model_list[0].genes]

        results = ind_data.copy()
        results['gene_ids'] = results['gene_ids'].apply(lambda x: list(set(x).intersection(genes)))
        results = results[results['gene_ids'].map(lambda x: len(x)) > 0]

        task_list = read_tasks('C:/Users/Sigve/Genome_Data/Human1/Human1_GEM/tasks/essential_tasks_min_ess_aa.tsv', model_list)
        # Adds empty task result list column
        results['tasks_results'] = np.empty((len(results), 0)).tolist()

        # Actual FBA
        results = parallelize_dataframe(results, partial(knockout_FBA_w_tasks, task_list, model_list), n_cores=14)

        results.loc[-1] = ['REF', [], round(model_list[0].slim_optimize(error_value='nan'), 3), [1]]

        # Put REF on top
        results.index = results.index + 1  # shifting index
        results.sort_index(inplace=True)

        results['gene_ids'] = results['gene_ids'].apply(';'.join)
        results['tasks_results'] = results['tasks_results'].apply(lambda x: x if not all(x) else ['ALL PASS'])

        results.reset_index(inplace=True, drop=True)
        results[['sample_ids', 'gene_ids', 'solution', 'tasks_results']].to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/results/ind_results/filtered/start_stop_het_inc_ess/ind_{0}.tsv'.format(tissue), sep='\t')

    end_time = time.time()

    print('Total time: %.6f seconds' % (end_time - start_time))


if __name__ == '__main__':
    main()
