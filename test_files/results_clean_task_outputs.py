
import pandas as pd

# make clean outputs. Flips results sideways so it is easy to see which tasks failed.

tissue_list = ['brain', 'pancreas', 'spleen', 'adipose_tissue', 'adrenal_gland', 'pituitary', 'thyroid', 'blood',
               'heart', 'kidney', 'liver', 'muscle', 'nerve', 'lung']
tissue_list.sort()

for tissue in tissue_list:

    input_file = 'C:/Users/Sigve/Genome_Data/results/ind_results/filtered/start_stop_het/full_tasks/ind_results_f_{0}.tsv'.format(tissue)
    output_file = 'C:/Users/Sigve/Genome_Data/results/ind_results/filtered/start_stop_het/full_tasks/ind_results_f_{0}_task_failures.tsv'.format(tissue)

    res = pd.read_table(input_file, index_col=0)
    res = res[res['tasks_results'] != "['ALL PASS']"]
    res['tasks_results'] = res['tasks_results'].apply(lambda x: [int(i) for i in x[1:-1].split(', ')])

    dict = {}
    for i, data in res.iterrows():
        dict[data['sample_ids']] = data['tasks_results']

    tasks_df_r = pd.DataFrame(dict)
    #print(tasks_df_r)

    file_path_t ='C:/Users/Sigve/Genome_Data/Human1/Human1_GEM/tasks/tissue_full_tasks/full_tasks_minus_ess_{0}.tsv'.format(tissue)
    tasks_df_t = pd.read_table(file_path_t, index_col=0).reset_index()
    tasks_df_t['task_id'] = tasks_df_t['index']#
    tasks_df_t['_'*89 + 'description'] = tasks_df_t['description'].apply(lambda x:  '_'*(100-len(x)) + x)
    tasks_df_t = tasks_df_t[['task_id', '_'*89 + 'description']]

    tasks_df_t = tasks_df_t.merge(tasks_df_r, right_index=True, left_index=True)

    tasks_df_t.to_csv(path_or_buf=output_file, sep='\t')

