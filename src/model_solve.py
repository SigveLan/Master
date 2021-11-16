
import re
import pandas as pd
import itertools
import cobra.flux_analysis
from multiprocessing import Pool
from functools import partial
import time
from assorted_functions import combinations_subset, parallelize_dataframe, knockout_FBA

"""A mess of a document with different code cells.
Good to to use for any testing that involves the Recon3D model as it takes some time to load in."""


def read_combinations(combinations_path: str) -> pd.DataFrame:
    """Prepare SNP combinations, reaction constraints"""

    combinations = pd.read_table(combinations_path, index_col=0)
    combinations['combinations'] = combinations['combinations'].apply(lambda x: x.split(';'))

    combinations['gene_model_ids'] = combinations['combinations']\
        .apply(lambda x: SNP_results.loc[SNP_results['variant_name'].isin(x), ['gene_number']].iloc[:, 0].tolist())

    id_list = ' '.join(model.genes.list_attr('id'))

    combinations['gene_model_ids'] = combinations['gene_model_ids'].apply(lambda x:
                                    list(set(itertools.chain.from_iterable(
                                    [re.findall(r"(?:\s)(" + str(i) + r"\w\w\w\d)", id_list) for i in x]))))

    return combinations


if __name__ == '__main__':
    start_time = time.time()
    model = cobra.io.load_json_model('C:/Users/Sigve/Genome_Data/Recon3D/JSON/Recon3D.json')
    SNP_results = pd.read_table('C:/Users/Sigve/Genome_Data/results/SNPs_effect_mod.tsv', index_col=0)
    end_time = time.time()
    print('Model load time: %.6f seconds' % (end_time - start_time))

    combinations = read_combinations('C:/Users/Sigve/Genome_Data/results/SNP_combinations.tsv')

    start_time = time.time()
    combinations = parallelize_dataframe(combinations, partial(combinations_subset, partial(knockout_FBA, model)))

    end_time = time.time()
    print('FBA run time: %.6f seconds' % (end_time - start_time))

    combinations.to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/results/test.tsv', sep='\t')
