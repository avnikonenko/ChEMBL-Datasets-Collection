import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import pandas as pd
import Dataset_collection as Dataset10
from multiprocessing import Pool

from functools import partial


def cuttof_target_type(targ):
    targ = targ.replace(' ','_').lower()

    if any(i in targ for i in ['epigenetic_factor','transporter']):
        ## plog_value
        act_dict = {'Ki': (6, '>='), 'IC50': (6, '>='), 'Kd': (6, '>='), 'EC50': (6, '>=')}
        inact_dict = {'Ki': (6, '<'), 'IC50': (6, '<'), 'Kd': (6, '<'), 'EC50': (6, '<')}

        type_act = ['allosteric', 'inverse agonist','antagonist', 'agonist', 'inhibitor', 'binding']
        contr_list = [('agonist', 'antagonist')]

    elif 'ion_channel' in targ:
        type_act = ['allosteric', 'inverse agonist', 'antagonist',
                    'agonist', 'binding', 'blockator']
        contr_list = [('agonist', 'antagonist')]
        # plog_value
        act_dict = {'Ki': (5, '>='), 'IC50': (5, '>='), 'Kd': (5, '>='), 'EC50': (5, '>=')}
        inact_dict = {'Ki': (5, '<'), 'IC50': (5, '<'), 'Kd': (5, '<'), 'EC50': (5, '<')}


    elif 'enzyme' in targ:
        type_act = ['allosteric', 'activator', 'inhibitor', 'binding']
        contr_list = [('activator', 'inhibitor')]
        # plog_value
        act_dict = {'Ki': (7.5, '>='), 'IC50': (7.5, '>='), 'Kd': (7.5, '>='), 'EC50': (7.5, '>=')}
        inact_dict = {'Ki': (7.5, '<'), 'IC50': (7.5, '<'), 'Kd': (7.5, '<'), 'EC50': (7.5, '<')}


    elif 'membrane_protein' or 'membrane_receptor' in targ:
        type_act = ['allosteric', 'inverse agonist','antagonist','agonist', 'inhibitor', 'binding']
        contr_list = [('agonist', 'antagonist'), ('inhibitor', 'agonist')]
        # plog_value
        act_dict = {'Ki': (7, '>='), 'IC50': (7, '>='), 'Kd': (7, '>='), 'EC50': (7, '>=')}
        inact_dict = {'Ki': (7, '<'), 'IC50': (7, '<'), 'Kd': (7, '<'), 'EC50': (7, '<')}

    else:
        # plog_value
        act_dict = {'Ki': (6, '>='), 'IC50': (6, '>='), 'Kd': (6, '>='), 'EC50': (6, '>=')}
        inact_dict = {'Ki': (6, '<'), 'IC50': (6, '<'), 'Kd': (6, '<'), 'EC50': (6, '<')}

        type_act = ['allosteric', 'inverse agonist','antagonist','agonist', 'inhibitor', 'binding']
        contr_list = [('agonist', 'antagonist'), ('inhibitor', 'agonist')]


    return act_dict, inact_dict, type_act, contr_list


def pool_map_iter(chembl_row, inp_path, type_dataset, input_fname_target, db, smi_std, sep):
    act_dict, inact_dict, type_act, contr_list, exp_type = None, None, None, None, None

    if 'target_class' in chembl_row:
        act_dict, inact_dict, type_act, contr_list = cuttof_target_type(chembl_row['target_class'])

    if 'cutoff_act' in chembl_row:
        act_dict = {i: (float(chembl_row['cutoff_act']), '>=') for i in chembl_row['activity'].split(';')},
        inact_dict = {i: (float(chembl_row['cutoff_inact']), '<') for i in chembl_row['activity'].split(';')},

    if 'category' in chembl_row:
        type_act = chembl_row['category'].split(';')
        contr_list = [('agonist', 'antagonist'), ('inhibitor', 'agonist'), ('activator', 'inhibitor'),
                      ('activator', 'antagonist')]

    if type_dataset == 'class':
        value = 'plog_value'
    else:
        value = 'value'
        exp_type = ['Ki', 'IC50', 'Kd', 'EC50']

    outfname = os.path.join(inp_path, type_dataset, value, input_fname_target, chembl_row['chembl_id'], chembl_row['chembl_id'])

    Dataset10.start(Chemblid=chembl_row['chembl_id'],
                    outfname=outfname,
                    type_act=type_act,
                    contr_list=contr_list,
                    act_dict=act_dict,
                    inact_dict=inact_dict,
                    value=value, db=db, smi_std=smi_std,
                    type_dataset=type_dataset,
                    exp_type=exp_type, sep=sep)


def get_from_csv(fname, type_dataset='class', sep='\t', ncpu=1, db=None, smi_std=None):

    data = pd.read_csv(fname, sep=sep)
    inp_path = os.path.dirname(fname)
    input_fname_target = os.path.basename(fname).split('.')[0]
    p = Pool(ncpu)

    p.map(partial(pool_map_iter,inp_path=inp_path, type_dataset=type_dataset,
                  input_fname_target=input_fname_target, db=db, smi_std=smi_std, sep=sep), [i[1] for i in data.iterrows()])



if __name__ == '__main__':
    parser = ArgumentParser(description='', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='target_class_file.txt', required=True)
    parser.add_argument('-t', '--type_dataset', metavar='class/reg', default='class')
    parser.add_argument('-d', '--db', metavar='ChEMBL-db.db', default=None)
    parser.add_argument('-s', '--smi_std', metavar='smi_std.smi', default=None)
    parser.add_argument('-n', '--ncpu', metavar='int', default=1, type=int)
    parser.add_argument('--sep', metavar=',', default='\t')

    args = parser.parse_args()
    print(vars(args))

    get_from_csv(fname=args.input, type_dataset=args.type_dataset, ncpu=args.ncpu, db=args.db, smi_std=args.smi_std, sep=args.sep)

