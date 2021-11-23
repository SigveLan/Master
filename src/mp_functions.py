import pandas as pd
import numpy as np
import cobra.flux_analysis
from multiprocessing import Pool, cpu_count

# File that contains functions pertaining to multiprocessing.

def split_filter(SNP_filter, SNPs_df: pd.DataFrame) -> pd.DataFrame:
    """Applies the SNP_filter function to a given subsets of the SNP data"""
    return SNP_filter(SNPs_df)


def knockout_FBA(model: cobra.Model, gene_ids: list) -> cobra.Solution:
    """Knock out FBA of given combination of genes."""
    with model:
        for gene_id in gene_ids:
            model.genes.get_by_id(gene_id).knock_out()
        return model.optimize()


def combinations_subset(knockout_FBA, combinations: pd.DataFrame) -> pd.DataFrame:
    """Applies knockout_FBA to the given data subset."""
    combinations['results'] = combinations['gene_model_ids'].apply(knockout_FBA)
    return combinations


def parallelize_dataframe(df: pd.DataFrame, func, n_cores=cpu_count()) -> pd.DataFrame:
    """Splits a dataframe into subsets and applies the given function, dividing the load over multiple CPU threads."""
    df_split = np.array_split(df, n_cores)
    pool = Pool(n_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df

