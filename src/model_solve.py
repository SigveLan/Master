
import re
import pandas as pd
import itertools
import cobra.flux_analysis
from functools import partial
import time
from mp_functions import combinations_subset, parallelize_dataframe, knockout_FBA

"""A file intended to be used for the model solving part instead of the jupyter notebook in the end."""


def apply_gene_ids_to_combinations(model: cobra.Model, SNP_results: pd.DataFrame, combinations: pd.DataFrame) -> pd.DataFrame:
    """Prepare SNP combinations for use as reaction constraints by producing model gene id lists for each combination"""

    combinations['combinations'] = combinations['combinations'].apply(lambda x: x.split(';'))

    combinations['gene_model_ids'] = combinations['combinations']\
        .apply(lambda x: SNP_results.loc[SNP_results['variant_name'].isin(x), ['gene_number']].iloc[:, 0].tolist())

    id_list = ';' + ';'.join(model.genes.list_attr('id'))

    combinations['gene_model_ids'] = combinations['gene_model_ids'].apply(lambda x:
                                        list(set(itertools.chain.from_iterable(
                                        [re.findall(r"(?:;)(" + str(i) + r"_AT\d)", id_list) for i in x]))))

    return combinations


def model_solve(SNPs_mod: pd.DataFrame, model_path: str, combinations_path: str, n_cores: int) -> pd.DataFrame:

    start_time = time.time()
    model = cobra.io.load_json_model(model_path)
    end_time1 = time.time()
    print('Model load time: %.6f seconds' % (end_time1 - start_time))

    combinations = pd.read_table(combinations_path, index_col=0)
    combinations = apply_gene_ids_to_combinations(model, SNPs_mod, combinations)

    print("KNock out FBA runs on " + str(combinations.shape[0]) + " gene combinations, divided over " +
          str(n_cores) + "CPU threads.")

    combinations = parallelize_dataframe(combinations, partial(combinations_subset, partial(knockout_FBA, model)), n_cores)
    end_time2 = time.time()
    print('FBA run time: %.6f seconds' % (end_time2 - end_time1))

    return combinations


