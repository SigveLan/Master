import itertools
import pandas as pd
import cobra.flux_analysis
from cobra import Metabolite, Reaction, Model
import time
import numpy as np
from functools import partial
from src.mp_functions import combinations_subset, parallelize_dataframe, knockout_FBA, knockout_FBA_w_tasks

from functools import partial
from src.met_task_functions import get_met_ids, constrain_model, create_reactions


def main ():


    start_time = time.time()
    model_file_path = 'C:/Users/Sigve/Genome_Data/Human1/Human1_GEM/GTEx/brain.xml'
    model = cobra.io.read_sbml_model(model_file_path)

    end_time = time.time()
    print('Model load time: %.6f seconds' % (end_time - start_time))

    model_list = constrain_model(model, ALLMETSIN=True)

    tasks_df = pd.read_table('C:/Users/Sigve/Genome_Data/Human1/Human1_GEM/tasks/essential_tasks.tsv')

    # Formatting data
    for b in ['LBin', 'LBout', 'UBin', 'UBout']:
        tasks_df[b] = tasks_df[b].apply(lambda x: x.split(','))

    for put in ['inputs', 'outputs']:
        tasks_df[put] = tasks_df[put].apply(lambda x: [e + ']' for e in x[1:-1].split(']')][0:-1])

    tasks_df['equations'] = tasks_df['equations'].apply(str)

    tasks_df[['met_ids', 'model_num']] = tasks_df.apply(partial(get_met_ids, model_list), axis=1, result_type='expand')

    tasks_df = create_reactions(tasks_df)




if __name__ == '__main__':
    main()
