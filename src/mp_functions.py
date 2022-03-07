import pandas as pd
import numpy as np
import cobra.flux_analysis
from cobra import Metabolite
from multiprocessing import Pool, cpu_count

# File that contains functions used for multiprocessing.


def split_filter(SNP_filter, SNPs_df: pd.DataFrame) -> pd.DataFrame:
    """Applies the SNP_filter function to a given subsets of the SNP data"""
    return SNP_filter(SNPs_df)


def knockout_FBA(model: cobra.Model, gene_ids: list) -> cobra.Solution:
    """Knock out FBA of given combination of genes."""
    with model:
        for gene_id in gene_ids:
            try:
                model.genes.get_by_id(gene_id).knock_out()
            except KeyError:
                return gene_id + ' not in model.'
        return model.optimize()


def knockout_FBA_w_tasks(tasks_df: pd.DataFrame, model_list: list, gene_ids: list) -> list:
    """Performs knockout FBA and checks tasks for the knockout."""
    with model_list[0]:
        for gene_id in gene_ids:
            try:
                model_list[0].genes.get_by_id(gene_id).knock_out()
            except KeyError:
                return gene_id + ' not in model.'
        res = [model_list[0].optimize().objective_value]

    for ind, data in tasks_df.iterrows():
        t_model = model_list[data.model_num]

        with t_model:
            for subset in [data.in_rx, data.out_rx]:
                for rx in subset:
                    if rx == 'ALLMETSIN':
                        # Adds boundary metabolites for other reactions when ALLMETSIN is used
                        for r in subset[1:]:
                            for m2 in r.metabolites:
                                for r2 in m2.reactions:
                                    if r2.boundary and r2.id != r.id:
                                        r2.add_metabolites({Metabolite(
                                                            m2.id[:-4] + 'x[x]',
                                                            formula=m2.formula,
                                                            name=' '.join(m2.name.split(' ')[:-1]) + ' [Boundary]',
                                                            compartment='x'): 1})
                        continue
                    t_model.add_reaction(rx)

            if data.equ != 'nan':
                t_model.add_reaction(data.equ)

            for gene_id in gene_ids:
                t_model.genes.get_by_id(gene_id).knock_out()

            if t_model.optimize().objective_value is None:
                res.append(0)
            else:
                res.append(1)
    return res


def combinations_subset(knockout_FBA, combinations: pd.DataFrame) -> pd.DataFrame:
    """Applies knockout_FBA to the given data subset."""
    combinations['results'] = combinations['gene_ids'].apply(knockout_FBA)
    return combinations


def parallelize_dataframe(df: pd.DataFrame, func, n_cores=cpu_count()) -> pd.DataFrame:
    """Splits a dataframe into subsets and applies the given function, dividing the load over multiple CPU threads."""
    df_split = np.array_split(df, n_cores)
    pool = Pool(n_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df

