
import itertools
import pandas as pd
import cobra.flux_analysis
from cobra import Metabolite
import time
import numpy as np
from src.mp_functions import combinations_subset, parallelize_dataframe, knockout_FBA, knockout_FBA_w_tasks, knockout_FBA_simple

from src.FBA_scripts.met_task_functions import constrain_model, read_tasks


tissues = ['spleen', 'adipose_tissue', 'adrenal_gland', 'pituitary', 'thyroid', 'kidney', 'liver', 'nerve', 'lung', 'skin']

for tissue in tissues:
    model_file_path = 'C:/Users/Sigve/Genome_Data/Human1/Human1_GEM/GTEx/{0}.xml'.format(tissue)
    model_list = constrain_model(model_file_path, ALLMETSIN_OUT=False)
    task_list = read_tasks('C:/Users/Sigve/Genome_Data/Human1/Human1_GEM/tasks/full_tasks_minus_ess.tsv', model_list)

    initial_res = pd.DataFrame([['REF', [], 0, []]], columns=['sample_ids', 'gene_ids', 'solution', 'tasks_results'])
    initial_res = knockout_FBA_w_tasks(task_list, model_list, initial_res)

    successful_tasks = initial_res.iat[0, 3]
    task_list = [task for task, p in zip(task_list, successful_tasks) if bool(p)]

    file_path_t = 'C:/Users/Sigve/Genome_Data/Human1/Human1_GEM/tasks/full_tasks_minus_ess.tsv'
    tasks_df_t = pd.read_table(file_path_t, index_col=0)


    tasks_df_t['REF'] = successful_tasks
    #tasks_df_t = tasks_df_t[tasks_df_t['REF'] == 1][['REF', 'description']]
    tasks_df_t.to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/Human1/Human1_GEM/tasks/tissue_full_tasks/full_tasks_minus_ess_{0}.tsv'.format(tissue), sep='\t')
