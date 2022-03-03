import cobra
import pandas as pd
from cobra import Metabolite, Reaction


def produce_task_models(task_file_path: str, model: cobra.Model) -> pd.DataFrame():
    """Creates one model for all tasks, be weary of RAM usage."""
    def get_met_ids(task: pd.Series) -> dict:

        comps = {'s': 'Extracellular',
                 'p': 'Peroxisome',
                 'm': 'Mitochondria',
                 'c': 'Cytosol',
                 'l': 'Lysosome',
                 'r': 'Endoplasmic reticulum',
                 'g': 'Golgi apparatus',
                 'n': 'Nucleus',
                 'i': 'Inner mitochondria',
                 'x': 'Boundary'}

        names = ['inputs', 'outputs']
        met_list = [task['inputs'], task['outputs']]

        if task['equations'] != 'nan':
            met_list.append( [item for item in task['equations'].split(' ') if item[-1] == ']'])
            names.append('equ')

        mets = task['models'].metabolites

        met_ids = {}
        for name, sub_list in zip(names, met_list):
            for met in sub_list:

                if met[:-3] == 'ALLMETSIN':
                    if name == 'inputs':
                        met_ids[met] = 'ALLMETSIN_IN'
                    elif name == 'outputs':
                        met_ids[met] = 'ALLMETSIN_OUT'
                    else:
                        raise ValueError("'ALLMETSIN' cannot be in equation.")

                elif met[:-3] == 'ALLMETS':
                    met_ids[met] = 'ALLMETS'
                    raise ValueError("'ALLMETS' is not supported.")

                else:
                    comp = met[-2]
                    temp = met[:-3] + ' [{0}]'.format(comps[comp])
                    for m in mets:
                        if m.name == temp:
                            met_ids[met] = m
                            break
                    else:
                        # Failed to find
                        raise ValueError("Failed to find metabolite for met_name: " + temp)

        return met_ids

    def constrain_model(models_df: pd.DataFrame) -> None:

        for ind, data in models_df.iterrows():

            if any(v == 'ALLMETSIN_IN' for k, v in data.met_ids.items()):
                for rx in data.models.exchanges:
                    rx.lower_bound = -1000
                    rx.upper_bound = 0

            elif any(v == 'ALLMETSIN_OUT' for k, v in data.met_ids.items()):
                for rx in data.models.exchanges:
                    rx.lower_bound = 0
                    rx.upper_bound = 1000

            else:
                for rx in data.models.exchanges:
                    m = list(rx.metabolites.keys())[0]
                    boundary_met = Metabolite(m.id[:-4] + 'x[x]', formula=m.formula,
                                              name=' '.join(m.name.split(' ')[:-1]) + ' [Boundary]', compartment='x')
                    rx.add_metabolites({boundary_met: 1})

    def create_reactions(tasks: pd.DataFrame) -> pd.DataFrame:
        # Producing reactions based on tasks
        rx_list = []
        in_list = []
        out_list = []

        for ind, data in tasks.iterrows():

            for i, name, ml, lbs, ubs in zip([1, -1], ['in', 'out'], [data.inputs, data.outputs], [data.LBin, data.LBout],
                                             [data.UBin, data.UBout]):

                rxl = []
                for j, m, lb, ub in zip(range(len(ml)), ml, lbs, ubs):

                    if m[:9] == 'ALLMETSIN':
                        for m2 in ml[1:]:
                            for r in data['met_ids'][m2].reactions:
                                if r.boundary:
                                    r.add_metabolites({Metabolite(m2.id[:-4] + 'x[x]', formula=m2.formula,
                                        name=' '.join(m2.name.split(' ')[:-1]) + ' [Boundary]', compartment='x'): 1})
                        continue

                    rx = Reaction('ess_{0}_{1}_{2}'.format(ind+1, name, j))
                    rx.add_metabolites({data['met_ids'][m]: i})
                    rx.lower_bound = float(lb)
                    rx.upper_bound = float(ub)
                    rxl.append(rx)

                if name == 'in':
                    in_list.append(rxl)
                else:
                    out_list.append(rxl)

            if data.equations != 'nan':
                t = [[data['met_ids'][subsub] for subsub in sub.split(' ') if len(subsub) > 1] for sub in data.equations.split('=')]
                d = {}

                for i, ml in zip([-1, 1], t):
                    for m in ml:
                        d[m] = i

                rx = Reaction('ess_{0}'.format(ind + 1))
                rx.add_metabolites(d)
                rx.lower_bound = float(data.LBequ)
                rx.upper_bound = float(data.UBequ)
                rx.name = data.description

                rx_list.append(rx)

            else:
                rx_list.append('nan')

        return pd.DataFrame(list(zip(in_list, out_list, rx_list, tasks['models'].tolist())), columns=['in_rx', 'out_rx', 'equ', 'models'])

    def apply_rxns(models_df: pd.DataFrame) -> None:

        for ind, data in models_df.iterrows():

            for rx in data.in_rx + data.out_rx:
                data.models.add_reaction(rx)

            if data.equ != 'nan':
                data.models.add_reaction(data.equ)

    tasks_df = pd.read_table(task_file_path)

    # Formatting data
    for b in ['LBin', 'LBout', 'UBin', 'UBout']:
        tasks_df[b] = tasks_df[b].apply(lambda x: x.split(','))

    for put in ['inputs', 'outputs']:
        tasks_df[put] = tasks_df[put].apply(lambda x: [e + ']' for e in x[1:-1].split(']')][0:-1])

    tasks_df['equations'] = tasks_df['equations'].apply(str)

    tasks_df['models'] = tasks_df.apply(lambda x: model.copy(), axis=1)


    tasks_df['met_ids'] = tasks_df.apply(get_met_ids, axis=1)

    constrain_model(tasks_df[['models', 'met_ids']])
    tasks_df = create_reactions(tasks_df)
    apply_rxns(tasks_df)

    return tasks_df
