import cobra
import pandas as pd
from functools import partial
from cobra import Metabolite, Reaction


def get_met_ids(model_list: list, task: pd.Series) -> list:
    """Find model metabolite objects from task metabolites"""

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
        met_list.append([item for item in task['equations'].split(' ') if item[-1] == ']'])
        names.append('equ')

    model_num = 1

    if any(met[:-3] == 'ALLMETSIN' for met in met_list[0]):
        model_num = 2
    if any(met[:-3] == 'ALLMETSIN' for met in met_list[1]):
        if model_num == 2:
            raise ValueError("'ALLMETSIN' cannot be in both input and output.")
        model_num = 3

    met_ids = {}
    for name, sub_list in zip(names, met_list):
        for met in sub_list:

            if met[:-3] == 'ALLMETSIN':
                met_ids[met] = 'ALLMETSIN'

            elif met[:-3] == 'ALLMETS':
                met_ids[met] = 'ALLMETS'
                raise ValueError("'ALLMETS' is not supported.")

            else:
                comp = met[-2]
                temp = met[:-3] + ' [{0}]'.format(comps[comp])
                for m in model_list[model_num].metabolites:
                    if m.name == temp:
                        met_ids[met] = m
                        break
                else:
                    # Failed to find
                    raise ValueError("Failed to find metabolite for met_name: " + temp)

    return [met_ids, model_num]


def constrain_model(model: cobra.Model, ALLMETSIN:bool = False) -> list:

    bound_model = model.copy()
    bound_model.name = 'bound_model'

    for rx in bound_model.exchanges:
        m = list(rx.metabolites.keys())[0]
        boundary_met = Metabolite(m.id[:-4] + 'x[x]', formula=m.formula,
                                  name=' '.join(m.name.split(' ')[:-1]) + ' [Boundary]', compartment='x')
        rx.add_metabolites({boundary_met: 1})

    md_list = [model, bound_model]

    if ALLMETSIN:

        ub_model = model.copy()
        ub_model.name = "upper_bound_exchanges_model"
        for rx in ub_model.exchanges:
            rx.lower_bound = -1000
            rx.upper_bound = 0
        md_list.append(ub_model)

        lb_model = model.copy()
        lb_model.name = "lower_bound_exchanges_model"
        for rx in lb_model.exchanges:
            rx.lower_bound = 0
            rx.upper_bound = 1000
        md_list.append(lb_model)

    return md_list


def create_reactions(tasks: pd.DataFrame) -> pd.DataFrame:
    """Producing reactions based on tasks."""

    rx_list = []
    in_list = []
    out_list = []

    for ind, data in tasks.iterrows():

        for i, name, ml, lbs, ubs in zip([1, -1], ['in', 'out'], [data.inputs, data.outputs], [data.LBin, data.LBout],
                                         [data.UBin, data.UBout]):

            rxl = []
            for j, m, lb, ub in zip(range(len(ml)), ml, lbs, ubs):

                if m[:9] == 'ALLMETSIN':
                    rxl.append('ALLMETSIN')
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

            for i, met_list in zip([-1, 1], t):
                for met in met_list:
                    d[met] = i

            rx = Reaction('ess_{0}'.format(ind + 1))
            rx.add_metabolites(d)
            rx.lower_bound = float(data.LBequ)
            rx.upper_bound = float(data.UBequ)
            rx.name = data.description

            rx_list.append(rx)

        else:
            rx_list.append('nan')

    return pd.DataFrame(list(zip(in_list, out_list, rx_list, tasks['model_num'].tolist())),
                        columns=['in_rx', 'out_rx', 'equ', 'model_num'])


def read_tasks(file_path: str, model_list: list) -> list:
    """Reads in metabolic tasks from a file, make sure the entries are in the correct form.
    Returns a list of lists containing reactions to be added to the model for each task."""

    tasks_df = pd.read_table(file_path)

    for b in ['LBin', 'LBout', 'UBin', 'UBout']:
        tasks_df[b] = tasks_df[b].apply(lambda x: x.split(','))

    for put in ['inputs', 'outputs']:
        tasks_df[put] = tasks_df[put].apply(lambda x: [e + ']' for e in x[1:-1].split(']')][0:-1])

    tasks_df['equations'] = tasks_df['equations'].apply(str)

    tasks_df[['met_ids', 'model_num']] = tasks_df.apply(partial(get_met_ids, model_list), axis=1, result_type='expand')

    tasks_df = create_reactions(tasks_df)

    return list(tasks_df.values.tolist())


