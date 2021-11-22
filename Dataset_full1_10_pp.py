from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import os.path
import sqlite3
from math import log10
from os import makedirs
from os.path import basename, exists, dirname

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from requests import get

url_chembl_assays = 'https://www.ebi.ac.uk/chembl/api/data/activity?format=json&limit=1000&offset={1}&target_chembl_id={0}'

"""
tree_act
class:self._type_act:keys[[any],[not],[all]]
"""
tree_act = {'allosteric': {'allosteric': [['allosteric'], [], []]},

            'antagonist': {'inverse agonist': [['inverse agonist'], [], []],
                           'competitive antagonist': [['competitive antagonist'], [], []],
                           'non-competitive antagonist': [
                               ['non-competitive', 'noncompetitive', 'non competitive', 'noncompetitively'],
                               ['inhibition'], []],
                           'antagonist': [['antagonist', 'antagonism', 'antagonistic activity', 'antagonistic'], [],
                                          []]},

            'agonist': {'partial agonist': [['partial agonist'], ['inhibition', 'inhibitory', 'inhibitor', 'antogonist',
                                                                  'inhibit', 'inactivation', 'intrinsic','inhibitors'], []],
                        'full agonist': [['maximal', 'maximum agonist'],
                                         ['inhibition', 'inhibitory', 'inhibitor','inhibitors', 'antogonist', 'inhibit',
                                          'inactivation', 'intrinsic'], ['maximum', 'stimulation']],

                        'agonist': [['agonist', 'agonistic activity', 'activity at', 'activation of',
                                     'activity towards', 'activity against', 'activity of', 'activation'],
                                    ['inhibition', 'inhibitory', 'inhibitor', 'inhibit', 'antagonist','inhibitors', 'inactivation'],
                                    []]},

            'activator': {'activator': [
                ['activation', 'activity at', 'activity towards', 'activity against', 'activity of', 'induction',
                 'activation'],
                ['inhibition', 'inhibitory','inhibitors', 'inhibitor', 'inhibit', 'antagonist', 'inactivation'], []]},

            'blockator': {'blockator': [['blockade', 'blockage'], [], []]},

            'inhibitor': {'competitive inhibitor': [['competitive'], [], []],
                          'uncompetitive inhibitor': [['uncompetitive'], [], []],
                          'non-competitive inhibitor': [['non-competitive', 'noncompetitive', 'non competitive'], [],
                                                        []],
                          'inhibitor': [['inhibitory activity', 'inhibitor', 'inhibitors', 'inhibition',
                                         'inhibitory concentration','activity of the inhibitors', 'inhibits', 'inhibit', 'inactivation'],
                                        [], []]},

            'binding': {'binding': [
                ['binding affinity', 'displacement', 'binding assay', 'competition', 'affinity', 'binding', 'displace'],
                ['inhibit', 'inhibition'], []]},
            'other': {'other': [[], [], []]}
            }


# list of search - right sequencing


# self._type_act = ['allosteric','competitive antagonist','inverse agonist',
#            'partial agonist', 'full agonist', 'antagonist', 'agonist',
#            'other']


class BaseDataSet:

    def __init__(self, filename=None, idchembl=None, type_act=None, overwrite=True, type_dataset='class'):
        """
        if filename:
        Warning: File must contains the next columns (necessary set):
         'cmp': ChEMBL ID;
         'smile':{smile structure};
         'operator': =, >=, <=, <, >, !=, ==;
         'value': value to assay;
         'units': unit to assay ;
         'bioactivity_type': type assay bioactivity;
         'pchembl_value': -log10(Value*10**-9)
         'act' : inhibitor, activator and etc.
         'type_act':
         'overwrite' - use only for given filedata. Rewrite default columns such as act, type act. True or False
        """
        if type_dataset not in ['class', 'reg']:
            print('Warning. Incorrect type dataset. Choose class or reg')

        self.__filename = filename
        self.__idchembl = idchembl

        self._id_cmp = 'cmp'

        self._type_act = type_act
        self._type_dataset = type_dataset
        self.overwrite = overwrite

        self._pdresult = pd.DataFrame()

        if self.__filename:
            self.__idchembl = None
            if not exists(self.__filename):
                raise ValueError("Error file")

    def load_data(self, local_db=None, sep='\t'):
        """
        local_db- path to the local SQLite chembl db (opt)
        """

        if self.__filename is not None and self.overwrite:
            res = pd.read_csv(self.__filename, sep=sep)
            res = self.__add_act_descr(res)
            self._pdresult = res
        else:
            if local_db is None:
                res = self.__webchembl()
            else:
                res = self.__load_local_chembl(local_db)
            if res is not None and not res.empty:
                res = self.__add_act_descr(res)
                self._pdresult = res
            else:
                print(res, self.__idchembl, 'Warning. Error loading')
                return False

        self._pdresult['target'] = self.__idchembl

        self._pdresult.loc[:,'operator'] = self._pdresult['operator'].fillna('=')
        self._pdresult.loc[:,'value'] = self._pdresult['value'].apply(pd.to_numeric, downcast='float', errors='coerce')
        self._pdresult.loc[:,'pchembl_value'] = self._pdresult['pchembl_value'].apply(pd.to_numeric, downcast='float',
                                                                                errors='coerce')
        self._pdresult = self._pdresult.dropna(subset=['value'])
        self._pdresult.loc[:, 'Log Note'] = None

        # convert uM or M to nM
        if 'uM' in self._pdresult['units'].tolist():
            print('uM')
            self._pdresult[self._pdresult['units'] == 'uM'] = self._pdresult[self._pdresult['units'] == 'uM'].apply(
                self.__convert_to_nmol, axis='columns')
        if 'M' in self._pdresult['units'].tolist():
            print('M')
            self._pdresult[self._pdresult['units'] == 'M'] = self._pdresult[self._pdresult['units'] == 'M'].apply(
                self.__convert_to_nmol, axis='columns')
        if True in self._pdresult['bioactivity_type'].isin(['pKd', 'pKi', 'pIC50', 'pEC50']):
            print('logValue')
            self._pdresult[self._pdresult['bioactivity_type'].isin(['pKd', 'pKi', 'pIC50', 'pEC50'])] = \
                self._pdresult[self._pdresult['bioactivity_type'].isin(['pKd', 'pKi', 'pIC50', 'pEC50'])].apply(
                    self.__convert_to_nmol, axis='columns')

        # pKi nM
        self._pdresult.loc[:,'plog_value'] = self._pdresult['value'].apply(
            lambda x: round(-log10(x * 10 ** -9), 2) if x > 0 else x)

        self._pdresult.loc[:,'pchembl_value'] = self._pdresult['pchembl_value'].round(2)

        self._pdresult.loc[:,'value'] = self._pdresult['value'].apply(lambda x: round(x, 3))

        if self._pdresult.empty:
            if self.__idchembl is not None:
                print(self.__idchembl, 'is not any results')
            else:
                print(self.__filename, 'is not any results')
            return False

        return True

    def __load_local_chembl(self, db_fname):
        """
        Description: Load and save all activity experiment result for target
        Arg: Target ChemblID
        Return: Compound set (description of act)

        """

        conn = sqlite3.connect(db_fname)
        cur = conn.cursor()
        cur.execute("""
            select 
              md.chembl_id, 
              act.standard_relation,
              act.standard_value,
              act.pchembl_value,
              act.standard_units, 
              act.standard_type, 
              ass.assay_type,
              ass.description,
              cs.canonical_smiles,
              d.journal || '-' || d.year || '-' || d.volume || '-' || d.issue || '-' || d.first_page || '-' || d.last_page || '. ' || d.doi
            from 
              target_dictionary as td, 
              assays as ass, 
              activities as act, 
              molecule_dictionary as md, 
              compound_structures as cs,
              docs as d
            where 
              td.tid = ass.tid and 
              ass.assay_id = act.assay_id and 
              act.molregno = md.molregno and 
              md.molregno = cs.molregno and
              ass.doc_id = d.doc_id and
              act.standard_value is not null and
              (ass.assay_type = 'B' or (act.standard_type in ('Ki','Kd')))  and
              ass.confidence_score >= 6 and
              act.data_validity_comment is null and
              td.chembl_id = ?
            ;
        """, (self.__idchembl,))
        res = cur.fetchall()
        res = pd.DataFrame(res, columns=['cmp', 'operator', 'value', 'pchembl_value', 'units', 'bioactivity_type',
                                         'assay_type',
                                         'assay_description', 'smile', 'reference'])
        return res

    def __webchembl(self):
        """

        Description: Load and save all activity experiment result for target
        Arg: Target ChemblID
        Return: Compound set (description of act)

        """
        # limit api - 1000 entr

        Next = True
        offset = 0
        cmps = []

        while Next:
            res = get(url_chembl_assays.format(self.__idchembl, offset))
            if not res.ok:
                print('Error conection')
                return None

            res = res.json()
            if res['page_meta']['next'] is None:
                Next = False
            offset += 1000

            for cmp in res['activities']:
                if cmp is None:
                    # or cmp['value'] is None or cmp['units'] is None or 'assay_description' not in cmp:
                    # print(cmp)
                    continue
                cmps.append({'cmp': cmp['molecule_chembl_id'],
                             'molecule_pref_name': cmp['molecule_pref_name'],
                             'operator': cmp['relation'],
                             'value': cmp.get('standart_value', cmp['value']),
                             'pchembl_value': cmp['pchembl_value'],
                             'units': cmp['units'],
                             'bioactivity_type': cmp['type'],
                             'assay_type': cmp['assay_type'],
                             'target_chembl_id': cmp['target_chembl_id'],
                             'target_pref_name': cmp['target_pref_name'],
                             'assay_description': cmp['assay_description'],
                             'smile': cmp['canonical_smiles'],
                             'reference': '-'.join([cmp['target_organism'], str(cmp['document_year']),
                                                    str(cmp['document_chembl_id'])])})
        pdresult = pd.DataFrame.from_dict(cmps)

        return pdresult

    def __getkeysearch(self, act_value):
        """

        return keys for search for type of act

        """
        act, keys = None, None
        for k, v in tree_act.items():
            if act_value in v:
                act = k
                keys = v[act_value]
                return act, keys

        return act, keys

    def __add_act_descr(self, res):

        res.loc[:,'act'] = 'other'
        res.loc[:,'type_act'] = 'other'

        for act_value in self._type_act:
            act, keys = self.__getkeysearch(act_value)
            if act is None and keys is None:
                print('Warning. Error act_value', act_value)
                continue

            # check if overwrite act
            m2 = res['act'] == 'other'
            # check assay
            m1 = res['assay_description'].apply(self.__check_act, args=(keys,))
            m = m1 & m2

            res.loc[m, 'act'] = act
            res.loc[m, 'type_act'] = act_value

        return res

    def __check_act(self, df, keys):
        line = ' {} '.format(df.lower()).replace('.', ' ').replace(',', ' ')
        if (any(' {} '.format(key) in line for key in keys[0]) or \
            (all(' {} '.format(key) in line for key in keys[2]) and keys[2])) and \
                not any(' {} '.format(key) in line for key in keys[1]):
            return True

        return False

    def add_standardized_smiles(self, smi_filename, colname='standardized_canonical_smiles'):
        if colname in self._pdresult.columns:
            print('Error colname {} already exist'.format(colname))
            if colname == 'standardized_canonical_smiles':
                self._pdresult = self._pdresult.drop(columns=[colname])
                print('Default colname. Continue')

            else:
                return None

        self._id_cmp = colname
        standart_smile = pd.read_csv(smi_filename, sep='\t', header=None, names=[colname, 'cmp'])
        self._pdresult = pd.merge(self._pdresult, standart_smile, on='cmp')

    def __convert_to_nmol(self, pdata):
        # error pandas/ twise apply for the first object
        if pdata['units'] == 'uM':
            coef = 10 ** 3
        elif pdata['units'] == 'M':
            coef = 10 ** 9

        elif pdata['bioactivity_type'] in ['pKd', 'pIC50', 'pEC50', 'pKi']:
            pdata['Log Note'] = ' '.join(['Original:', str(pdata['value']), pdata['bioactivity_type']])
            pdata['units'] = 'nM'
            pdata['bioactivity_type'] = pdata['bioactivity_type'].replace('p', '')
            pdata['value'] = 10 ** -pdata['value'] * 10 ** 9
            return pdata
        else:
            return pdata

        pdata['Log Note'] = ' '.join(['Original:', str(pdata['value']), pdata['units']])
        pdata['units'] = 'nM'
        # pdata['pchembl_value'] = round(-log10(pdata['value'] * 10 ** -6), 1) if pdata['value'] > 0 else pdata['value']
        pdata['value'] = pdata['value'] * coef

        return pdata

    def calc_MW(self):
        def _calc_MW(smi):
            mol = Chem.MolFromSmiles(smi)
            return Descriptors.MolWt(mol)

        self._pdresult.loc[:, 'MW'] = self._pdresult['standardized_canonical_smiles'].apply(_calc_MW)

    def save_data(self, arg='log', outfile=None, na_filter=True, drop_duplicate=True, make_dir=False, sep='\t'):
        """
        Input arg:

         - 'log'  - save log file for all results and don't filter status
         - 'act' - split of act for compound (filter for undefined status)
        (agonist, antagonist, allosteric, other)
         - 'type_act' - split of type act for compound (filter for undefined status)
        ( 'allosteric','competitive antagonist','inverse agonist','partial agonist', 'full agonist', 'antagonist', 'agonist', 'other')
        - 'all' - save all variant

        -outfile - preffics for result file. Default - '{filename}_{key_data}.csv'

        """
        if arg not in ['log', 'all', 'act', 'type_act']:
            print('Incorrect arg. Choose log/all/act/type_act')
            return False

        if not outfile:
            if self.__filename:
                f = basename(self.__filename).rsplit('.', 1)[0]
                work_path = dirname(self.__filename)
            else:
                f = self.__idchembl
                work_path = os.getcwd()
            outfile = '{file}'.format(file=f)
        else:
            work_path = dirname(outfile)
            outfile = basename(outfile).rsplit('.', 1)[0]
            if work_path and make_dir:
                makedirs(work_path, exist_ok=True)

        print('Target: ', outfile, arg)

        if arg == 'log' or arg == 'all':
            self._pdresult.to_csv(os.path.join(work_path, '{0}_{1}_log'.format(outfile, self._type_dataset)),
                                  sep=sep, index=False)

        if arg == 'act' or arg == 'all':
            self.__save_split_group('act', tree_act, os.path.join(work_path, outfile), na_filter, drop_duplicate,
                                    sep=sep)

        if arg == 'type_act' or arg == 'all':
            self.__save_split_group('type_act', self._type_act, os.path.join(work_path, outfile), na_filter,
                                    drop_duplicate, sep=sep)

    def __save_split_group(self, col, val, outfile, na_filter=True, drop_duplicate=True, sep='\t'):
        stat = pd.DataFrame()
        for key in val:
            res = pd.DataFrame()

            if self._type_dataset == 'class':
                res = self._pdresult[(self._pdresult[col] == key) &  # (self._pdresult.status != 'undefined') &
                                     (self._pdresult.status != 'controversial') & (self._pdresult.status != 'del') & (
                                             self._pdresult['exp_type'] == 'norm')]

            elif self._type_dataset == 'reg':
                res = self._pdresult[(self._pdresult[col] == key) & (self._pdresult['diff'] == 'norm') & (
                        self._pdresult['exp_type'] == 'norm')]

            if res.empty:
                continue

            if drop_duplicate:
                res = res.drop_duplicates(self._id_cmp)
                res = res.drop_duplicates('cmp')

            if na_filter:
                res.dropna(subset=[self._id_cmp]).to_csv('{0}_{1}_{2}.csv'.format(outfile, key, self._type_dataset),
                                                         sep=sep, index=False)
            else:
                res.to_csv('{0}_{1}_{2}.csv'.format(outfile, key, self._type_dataset), sep=sep, index=False)

            if self._type_dataset == 'class':
                stat.loc[key,'active'] = len(res[res['status'] == 'active'])
                stat.loc[key, 'inactive'] = len(res[res['status'] == 'inactive'])
                stat.loc[key, 'undefined'] = len(res[res['status'] == 'undefined'])

        if not stat.empty:
            stat = stat.rename_axis(index='act').reset_index().astype('int')
            stat.to_csv('{0}_{1}'.format(outfile, 'stat'), sep=sep, index=False)


class ReggDataSet(BaseDataSet):
    def __init__(self, filename=None, idchembl=None,
                 type_act=['allosteric', 'competitive antagonist', 'inverse agonist',
                           'partial agonist', 'full agonist', 'antagonist', 'agonist',
                           'inhibitor', 'other']):
        super().__init__(filename=filename, idchembl=idchembl, type_act=type_act, type_dataset='reg')

    def filt_data(self, bioactivity_type=['Ki', 'IC50', 'Kd'], units=['nM'], operator=['='], delta_diff=0.5):
        self._filt_experiment_type(bioactivity_type=bioactivity_type, units=units, operator=operator)
        self._filt_diff_results(delta_diff=delta_diff, operator=operator)

    def _filt_experiment_type(self, bioactivity_type=['Ki', 'IC50', 'Kd'], units=['nM'], operator=['=']):
        mask_true = self._pdresult['bioactivity_type'].isin(bioactivity_type) & self._pdresult['units'].isin(units) & \
                    self._pdresult['operator'].isin(operator)

        self._pdresult.loc[:, 'exp_type'] = 'undefined'
        self._pdresult.loc[mask_true, 'exp_type'] = 'norm'

    def _filt_diff_results(self, delta_diff=0.5, operator=['=']):
        def __mask_del_value(cmp_group, delta_diff):
            if len(cmp_group) > 1:
                if not all(cmp_group['operator'].isin(operator)) or \
                        cmp_group['plog_value'].max() - cmp_group['plog_value'].min() > delta_diff:
                    return True
            return False

        mask_del_diff = self._pdresult.groupby([self._id_cmp, 'bioactivity_type', 'units','act', 'type_act']).filter(
            __mask_del_value, delta_diff=delta_diff).index
        self._pdresult.loc[:, 'diff'] = 'norm'
        self._pdresult.loc[mask_del_diff, 'diff'] = 'del'

    def mean_duplicate(self):
        self._pdresult.loc[:, 'mean_plog_value'] = \
            self._pdresult.groupby([self._id_cmp, 'bioactivity_type','units', 'act', 'type_act'])[
                'plog_value'].transform(np.mean)
        self._pdresult.loc[:, 'mean_value'] = \
            self._pdresult.groupby([self._id_cmp, 'bioactivity_type', 'units', 'act', 'type_act'])[
                'value'].transform(np.mean)


class ClassDataSet(BaseDataSet):

    def __init__(self, filename=None, idchembl=None,
                 type_act=['allosteric', 'competitive antagonist', 'inverse agonist',
                           'partial agonist', 'full agonist', 'antagonist', 'agonist',
                           'inhibitor', 'other']):

        super().__init__(filename=filename, idchembl=idchembl, type_act=type_act, type_dataset='class')

    def add_status(self, filt={'Ki': (1000, '<='),
                               'IC50': (1000, '<='),
                               'Activity': (100, '>'),
                               'EC50': (1000, '<='),
                               'Kd': (1000, '<='),
                               'Inhibition': (70, '<='),
                               'FC': (0, '>'),
                               'Potency': (1000, '<=')},
                   status='active', unit=['%', 'nM'], value='value'):
        '''
        :param filt: dict {'Ki': (1000, '<='), 'IC50': (1000, '<=')} for value or {'KI':6} for plog_value
        :param status: active/inactive or some other
        :param unit: list ['nM','%']
        :param value: value (if filt value: Ki, IC50) or plog_value (if filt -log10(value): pKi, pIC50)
        :return:
        '''

        if 'status' not in self._pdresult.columns:
            self._pdresult['status'] = 'undefined'

        for bioactivity, cond in filt.items():
            val, match = cond
            mask = self.__get_index_act(bioactivity, val, match, status, unit, value)
            self._pdresult.loc[mask, 'status'] = status

    def __get_index_act(self, bioactivity, val, match, status, unit, value):

        mask1 = self._pdresult['bioactivity_type'] == bioactivity

        if match == '<=':
            mask2 = self._pdresult[value] <= float(val)
        elif match == '>=':
            mask2 = self._pdresult[value] >= float(val)
        elif match == '!=':
            mask2 = self._pdresult[value] != float(val)
        elif match == '==' or match == '=':
            mask2 = self._pdresult[value] == val
        elif match == '>':
            mask2 = self._pdresult[value] > float(val)
        else:  # '<'
            mask2 = self._pdresult[value] < float(val)

        mask3 = self._pdresult['units'].isin(unit)

        # nM, M, uM <- activity, -> inact
        if status == 'active':
            mask4 = self._pdresult['operator'].isin(['=', '<', '<=', '~'])
            return mask1 & mask2 & mask3 & mask4

        if status == 'inactive':
            mask4 = self._pdresult['operator'].isin(['=', '>', '>=', '~'])
            return mask1 & mask2 & mask3 & mask4

        return mask1 & mask2 & mask3

    def filt_data(self, bioactivity_type=['Ki', 'IC50', 'Kd'], units=['nM']):
        self._filt_experiment_type(bioactivity_type=bioactivity_type, units=units)
        self._filt_contr_status()

    def _filt_experiment_type(self, bioactivity_type=['Ki', 'IC50', 'Kd'], units=['nM']):
        mask_true = self._pdresult['bioactivity_type'].isin(bioactivity_type) & self._pdresult['units'].isin(units)

        self._pdresult.loc[:, 'exp_type'] = 'undefined'
        self._pdresult.loc[mask_true, 'exp_type'] = 'norm'

    def _filt_contr_status(self):
        """
        Check active and inactive status for one act
        """

        def __check_nonunique_status(df):
            stat = df['status'].unique()
            if df['act'].tolist()[0] != 'other':
                if len(stat) > 1:
                    # if 'active' in stat and 'inactive' in stat or undefined
                    return True
            return False

        mask_nonunique = self._pdresult.groupby([self._id_cmp, 'act']).filter(__check_nonunique_status).index
        self._pdresult.loc[mask_nonunique, 'status'] = 'del'

    def filt_contr_data(self, contr=tuple(('agonist', 'antagonist'))):
        self._filt_contr_act_result(contr)

    def _filt_contr_act_result(self, contr):
        def __check_contr_type_act(cmp_group, contr):
            if contr[0] in cmp_group['act'].unique() and contr[1] in cmp_group['act'].unique():
                return True
            return False

        contr_mask = self._pdresult[self._pdresult['status'] == 'active'].groupby(self._id_cmp).filter(
            __check_contr_type_act, contr=contr).index
        self._pdresult.loc[contr_mask, 'status'] = 'controversial'


def start(Chemblid, outfname, type_act, contr_list=None, act_dict=None, inact_dict=None, value='value', db=None, smi_std=None,
          type_dataset='class', filename=None, sep='\t', exp_type=None):
    '''
    :param Chemblid:
    :param outfname: fname
    :param type_act: list. Exp: ['other']
    :param contr_list:
    :param act_dict:
    :param inact_dict:
    :param value: value or plog_value
    :param db: chembldb fname
    :param smi_std: fname
    :param type_dataset: reg or class
    :param filename:
    :param sep:
    :param exp_type: list. Exp: ['Ki']
    :return:
    '''
    print('start', Chemblid, outfname, type_act, contr_list, act_dict, inact_dict, value, db, smi_std,
          type_dataset, filename, sep,exp_type)

    if type_dataset == 'class':
        if filename is not None:
            A = ClassDataSet(filename=filename, type_act=type_act)
        else:
            A = ClassDataSet(idchembl=Chemblid, type_act=type_act)
    else:
        if filename is not None:
            A = ReggDataSet(filename=filename, type_act=type_act)
        else:
            A = ReggDataSet(idchembl=Chemblid, type_act=type_act)

    if A.load_data(local_db=db, sep=sep):
        if smi_std is not None:
            A.add_standardized_smiles(smi_std, colname='standardized_canonical_smiles')
            print('id_cmp', A._id_cmp)
            A.calc_MW()

        if type_dataset == 'class':
            A.add_status(
                # act_dict = {'Ki': (100, '<='), 'IC50': (100, '<='), 'Kd': (100, '<=')},
                act_dict, 'active', unit=['nM'], value=value)
            A.add_status(
                # inact_dict = {'Ki': (10000, '>='), 'IC50': (10000, '>='), 'Kd': (10000, '>=')},
                inact_dict, 'inactive', unit=['nM'], value=value)
            A.add_status(
                # filt = {'Ki': (0, '='), 'IC50': (0, '='), 'Kd': (0, '=')},
                filt={k: (0, '=') for k in inact_dict.keys()},
                status='inactive', unit=['nM'], value=value)

            for contr in contr_list:
                A.filt_contr_data(contr)

            A.filt_data(bioactivity_type=list(act_dict.keys() | inact_dict.keys()), units=['nM'])
        else:
            A.mean_duplicate()
            A.filt_data(bioactivity_type=exp_type, units=['nM'], delta_diff=0.5)

        A.save_data(arg='type_act', outfile=outfname, make_dir=True, sep=sep)
        A.save_data(arg='log', outfile=outfname, make_dir=True, sep=sep)

        return True


if __name__ == '__main__':
    parser = ArgumentParser(description='', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='ChEMBLlID', required=True)
    parser.add_argument('--t', metavar='class/reg', required=True)
    parser.add_argument('--activity', metavar='Ki IC50 EC50', required=True, nargs='*')
    parser.add_argument('-o', '--output', metavar='output', required=False, default=None)
    parser.add_argument('-a', '--act', metavar='binding inhibitor activator', required=False, nargs='*', default=['other'])
    parser.add_argument('-c', '--contr', metavar='inhibitor activator', required=False, nargs='*', default=[])
    parser.add_argument('-v','--value', metavar='value/plog_value', required=False,  default='value')
    parser.add_argument('--active_value', metavar='7', required=False,  default='6')
    parser.add_argument('--active_op', metavar='>=', required=False,  default='>=')
    parser.add_argument('--inactive_value', metavar='7', required=False,  default='6')
    parser.add_argument('--inactive_op', metavar='<', required=False,  default='<')
    parser.add_argument('--sep', metavar=',/;', required=False,  default=';')

    args = parser.parse_args()
    print(vars(args))

    # for i in d['chemblid'].to_list():
    # start(i, 'Ki_1000_ready/'+i, ['other'],
    # contr_list=None, act_dict=None, inact_dict=None, value='plog_value',
    # db=chembl_23.db',
    # smi_std='chembl_23_step2_chembl_id.smi',
    # type_dataset='reg', sep=',', exp_type=['Ki'])

    _db = 'chembl_23.db'
    _smi_std = 'hembl_23_step2_chembl_id.smi'

    _chemblID = args.input
    _type_act = args.act
    _contr_list = args.contr
    _value = args.value
    _type_dataset = args.t
    _output = args.output

    _act_val = args.active_value
    _inact_val = args.inactive_value
    _act_op = args.active_op
    _inact_op = args.inactive_op

    _activity_list = args.activity
    _sep = args.sep

    if _output is None:
        _output = os.path.join(_chemblID,_chemblID)

    # _act_dict = {'Ki': (6, '>='), 'IC50': (6, '>='), 'Kd': (6, '>='), 'EC50': (6, '>=')}
    # _inact_dict = {'Ki': (6, '<'), 'IC50': (6, '<'), 'Kd': (6, '<'), 'EC50': (6, '<')}
    _act_dict = {i:(_act_val, _act_op) for i in _activity_list}
    _inact_dict = {i:(_inact_val, _inact_op) for i in _activity_list}


    start(Chemblid=_chemblID, outfname= _output, type_act=_type_act, value=_value,
          contr_list=_contr_list, act_dict=_act_dict, inact_dict=_inact_dict,
          db=_db, smi_std=_smi_std,type_dataset=_type_dataset,
          filename=None, sep=_sep, exp_type=_activity_list)


#cat chembl_id | xargs -I {1} Dataset_full1_10_pp.py -i {1} --t reg --activity EC50 IC50 -v plog_value -a other

    # salt problem -> smi_std

    # type act - warning order for search
    # type_act = ['allosteric', 'competitive antagonist', 'inverse agonist',
    #             'partial agonist', 'full agonist', 'antagonist', 'agonist',
    #             'inhibitor', 'other']
    # type_act = ['inhibitor', 'other']
    #  type_act = ['other']

    # use_target_fname('set.txt', db='../../chembl_23.db', type_dataset='class')
    # use_target_fname('set.txt', db='../../chembl_23.db', type_dataset='reg')
    # use_target_fname('set.txt', db='chembl_23.db', type_dataset='class',
    #                  smi_std='chembl_23_step2_chembl_id.smi')

# for i in d['chemblid'].to_list():
#     start(i, 'Ki_1000_ready/'+i, ['other'], contr_list=None, act_dict=None, inact_dict=None, value='plog_value', db='chembl_23.db',
#           smi_std='chembl_23_step2_chembl_id.smi',
#           type_dataset='reg', sep=',', exp_type=['Ki'])
