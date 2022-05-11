import pandas as pd
import numpy as np
import itertools
from multiprocessing import Pool, cpu_count


def combination_generation(temp_dfs: np.ndarray) -> pd.DataFrame:

    combinations_df = pd.DataFrame(columns=['phewas_code', 'gene_ids', 'snp_ids'])

    for temp_df in temp_dfs:
        code = temp_df['phewas_code'].iloc[0]
        SNP_ids = temp_df['variant_name'].tolist()
        combs = []
        # Set Range, could be input
        for L in range(5, 6 + 1):
            for subset in itertools.combinations(SNP_ids, L):
                combs.append(subset)
        gene_id_combs = []
        for subset in combs:
            gene_ids = []
            for i in subset:
                gene_ids.append(temp_df[temp_df['variant_name'] == i]['gene_id'].iloc[0].split('.')[0])
            gene_id_combs.append(';'.join(gene_ids))

        combinations_subset = pd.DataFrame({'snp_ids': combs, 'gene_ids': gene_id_combs})
        combinations_subset['snp_ids'] = combinations_subset['snp_ids'].apply(lambda x: ';'.join(x))

        combinations_subset = (combinations_subset.groupby('gene_ids').agg({'snp_ids': lambda x: x.tolist()})).reset_index()
        combinations_subset['phewas_code'] = code

        combinations_subset = combinations_subset[['phewas_code', 'gene_ids', 'snp_ids']]
        combinations_df = pd.concat([combinations_df, combinations_subset]).reset_index(drop=True)

    return combinations_df


def phewas_parallelization(df_list: list, func, n_cores=cpu_count()) -> pd.DataFrame:

    dfs = np.array_split(df_list, n_cores)
    pool = Pool(n_cores)
    df = pd.concat(pool.map(func, dfs)).reset_index(drop=True)
    pool.close()
    pool.join()
    return df
