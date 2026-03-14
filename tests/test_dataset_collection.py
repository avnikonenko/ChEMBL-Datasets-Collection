import pandas as pd
import numpy as np
import pytest

from Dataset_collection import ClassDataSet, ReggDataSet


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_class_ds(rows):
    """Return a ClassDataSet with _pdresult pre-populated from a list of dicts."""
    ds = ClassDataSet(type_act=['inhibitor', 'agonist', 'antagonist', 'other'])
    ds._pdresult = pd.DataFrame(rows)
    return ds


def _make_reg_ds(rows):
    ds = ReggDataSet(type_act=['inhibitor', 'other'])
    ds._pdresult = pd.DataFrame(rows)
    return ds


def _base_row(**kwargs):
    defaults = dict(
        cmp='CHEMBL1', smiles='C', operator='=', value=100.0,
        units='nM', bioactivity_type='Ki', pchembl_value=7.0,
        act='inhibitor', type_act='inhibitor',
        assay_description='inhibition assay', plog_value=7.0,
    )
    defaults.update(kwargs)
    return defaults


# ---------------------------------------------------------------------------
# __check_act  (keyword matching)
# ---------------------------------------------------------------------------

class TestCheckAct:
    """Tests for BaseDataSet.__check_act (accessed via name mangling)."""

    def setup_method(self):
        self.ds = ClassDataSet(type_act=['inhibitor'])
        self._check = self.ds._BaseDataSet__check_act

    def test_any_keyword_matches(self):
        keys = [['inhibition'], [], []]
        assert self.ds._BaseDataSet__check_act('Inhibition of kinase activity', keys) is True

    def test_exclusion_keyword_blocks_match(self):
        keys = [['agonist'], ['inhibition'], []]
        assert self._check('agonist inhibition assay', keys) is False

    def test_required_all_keywords_pass(self):
        # keys[2]: all must be present
        keys = [[], [], ['maximum', 'stimulation']]
        assert self._check('maximum stimulation of receptor', keys) is True

    def test_required_all_keywords_partial_fail(self):
        keys = [[], [], ['maximum', 'stimulation']]
        assert self._check('maximum activation of receptor', keys) is False

    def test_case_insensitive(self):
        keys = [['inhibitor'], [], []]
        assert self._check('INHIBITOR assay', keys) is True

    def test_punctuation_stripped(self):
        # dots and commas around keyword should still match
        keys = [['inhibitor'], [], []]
        assert self._check('an inhibitor. assay', keys) is True

    def test_empty_any_and_empty_all_returns_false(self):
        # keys[0] empty, keys[2] empty → no match condition satisfied
        keys = [[], [], []]
        assert self._check('any description here', keys) is False

    def test_partial_word_does_not_match(self):
        # 'inhibit' should not match keyword ' inhibitor ' due to space-padding
        keys = [['inhibitor'], [], []]
        assert self._check('inhibitory compound', keys) is False


# ---------------------------------------------------------------------------
# __convert_to_nmol  (unit conversion)
# ---------------------------------------------------------------------------

class TestConvertToNmol:
    def setup_method(self):
        self.ds = ClassDataSet(type_act=['inhibitor'])
        self._convert = self.ds._BaseDataSet__convert_to_nmol

    def _row(self, value, units, bioactivity_type='Ki'):
        return pd.Series(dict(value=value, units=units, bioactivity_type=bioactivity_type,
                              pchembl_value=None, act='inhibitor', type_act='inhibitor',
                              Log_Note=None))

    def test_uM_to_nM(self):
        row = self._row(1.0, 'uM')
        result = self._convert(row)
        assert result['value'] == pytest.approx(1000.0)
        assert result['units'] == 'nM'

    def test_M_to_nM(self):
        row = self._row(1e-6, 'M')
        result = self._convert(row)
        assert result['value'] == pytest.approx(1000.0)
        assert result['units'] == 'nM'

    def test_pKi_back_converts(self):
        row = self._row(7.0, 'nM', 'pKi')
        result = self._convert(row)
        assert result['value'] == pytest.approx(100.0, rel=1e-3)
        assert result['units'] == 'nM'
        assert result['bioactivity_type'] == 'Ki'

    def test_pIC50_back_converts(self):
        row = self._row(6.0, 'nM', 'pIC50')
        result = self._convert(row)
        assert result['value'] == pytest.approx(1000.0, rel=1e-3)
        assert result['bioactivity_type'] == 'IC50'

    def test_unknown_unit_unchanged(self):
        row = self._row(500.0, 'nM')
        result = self._convert(row)
        assert result['value'] == pytest.approx(500.0)
        assert result['units'] == 'nM'


# ---------------------------------------------------------------------------
# ClassDataSet.add_status
# ---------------------------------------------------------------------------

class TestAddStatus:
    def test_active_label(self):
        ds = _make_class_ds([_base_row(value=50.0, plog_value=7.3, units='nM', bioactivity_type='Ki', operator='=')])
        ds.add_status(filt={'Ki': (100.0, '<=')}, status='active', unit=['nM'], value='value')
        assert ds._pdresult.iloc[0]['status'] == 'active'

    def test_inactive_label(self):
        ds = _make_class_ds([_base_row(value=5000.0, plog_value=5.3, units='nM', bioactivity_type='Ki', operator='=')])
        ds.add_status(filt={'Ki': (1000.0, '>=')}, status='inactive', unit=['nM'], value='value')
        assert ds._pdresult.iloc[0]['status'] == 'inactive'

    def test_above_threshold_not_active(self):
        ds = _make_class_ds([_base_row(value=2000.0, units='nM', bioactivity_type='Ki', operator='=')])
        ds.add_status(filt={'Ki': (1000.0, '<=')}, status='active', unit=['nM'], value='value')
        assert ds._pdresult.iloc[0]['status'] == 'undefined'

    def test_wrong_unit_not_labeled(self):
        ds = _make_class_ds([_base_row(value=50.0, units='uM', bioactivity_type='Ki', operator='=')])
        ds.add_status(filt={'Ki': (100.0, '<=')}, status='active', unit=['nM'], value='value')
        assert ds._pdresult.iloc[0]['status'] == 'undefined'

    def test_active_requires_le_or_lt_operator(self):
        # operator '>' should not be labeled active
        ds = _make_class_ds([_base_row(value=50.0, units='nM', bioactivity_type='Ki', operator='>')])
        ds.add_status(filt={'Ki': (100.0, '<=')}, status='active', unit=['nM'], value='value')
        assert ds._pdresult.iloc[0]['status'] == 'undefined'

    def test_zero_value_labeled_inactive(self):
        ds = _make_class_ds([_base_row(value=0.0, units='nM', bioactivity_type='Ki', operator='=')])
        ds.add_status(filt={'Ki': (0, '=')}, status='inactive', unit=['nM'], value='value')
        assert ds._pdresult.iloc[0]['status'] == 'inactive'


# ---------------------------------------------------------------------------
# ClassDataSet._filt_contr_status  (contradictory active+inactive → del)
# ---------------------------------------------------------------------------

class TestFiltContrStatus:
    def _ds_with_status(self, rows):
        ds = _make_class_ds(rows)
        ds._pdresult['status'] = [r.get('status', 'undefined') for r in rows]
        ds._pdresult['exp_type'] = 'norm'
        return ds

    def test_active_and_inactive_same_compound_marked_del(self):
        rows = [
            _base_row(cmp='CHEMBL1', act='inhibitor', status='active'),
            _base_row(cmp='CHEMBL1', act='inhibitor', status='inactive'),
        ]
        ds = self._ds_with_status(rows)
        ds._filt_contr_status()
        assert all(ds._pdresult['status'] == 'del')

    def test_only_active_not_marked_del(self):
        rows = [_base_row(cmp='CHEMBL1', act='inhibitor', status='active')]
        ds = self._ds_with_status(rows)
        ds._filt_contr_status()
        assert ds._pdresult.iloc[0]['status'] == 'active'

    def test_other_act_not_marked_del(self):
        # act='other' compounds should be skipped even if active+inactive
        rows = [
            _base_row(cmp='CHEMBL1', act='other', status='active'),
            _base_row(cmp='CHEMBL1', act='other', status='inactive'),
        ]
        ds = self._ds_with_status(rows)
        ds._filt_contr_status()
        statuses = ds._pdresult['status'].tolist()
        assert 'del' not in statuses


# ---------------------------------------------------------------------------
# ClassDataSet._filt_contr_act_result  (controversial agonist+antagonist)
# ---------------------------------------------------------------------------

class TestFiltContrActResult:
    def test_active_in_both_contr_types_marked_controversial(self):
        rows = [
            _base_row(cmp='CHEMBL1', act='agonist', status='active'),
            _base_row(cmp='CHEMBL1', act='antagonist', status='active'),
        ]
        ds = _make_class_ds(rows)
        ds._pdresult['status'] = ['active', 'active']
        ds._pdresult['exp_type'] = 'norm'
        ds._filt_contr_act_result(('agonist', 'antagonist'))
        assert all(ds._pdresult['status'] == 'controversial')

    def test_active_in_only_one_type_unchanged(self):
        rows = [_base_row(cmp='CHEMBL1', act='agonist', status='active')]
        ds = _make_class_ds(rows)
        ds._pdresult['status'] = ['active']
        ds._pdresult['exp_type'] = 'norm'
        ds._filt_contr_act_result(('agonist', 'antagonist'))
        assert ds._pdresult.iloc[0]['status'] == 'active'


# ---------------------------------------------------------------------------
# ReggDataSet._filt_diff_results
# ---------------------------------------------------------------------------

class TestFiltDiffResults:
    def _reg_ds(self, rows):
        ds = _make_reg_ds(rows)
        ds._pdresult['exp_type'] = 'norm'
        ds._pdresult['diff'] = 'norm'
        return ds

    def test_large_range_marked_del(self):
        rows = [
            _base_row(cmp='CHEMBL1', act='inhibitor', type_act='inhibitor',
                      plog_value=7.0, operator='=', units='nM', bioactivity_type='Ki'),
            _base_row(cmp='CHEMBL1', act='inhibitor', type_act='inhibitor',
                      plog_value=7.8, operator='=', units='nM', bioactivity_type='Ki'),
        ]
        ds = self._reg_ds(rows)
        ds._filt_diff_results(delta_diff=0.5)
        assert all(ds._pdresult['diff'] == 'del')

    def test_small_range_stays_norm(self):
        rows = [
            _base_row(cmp='CHEMBL1', act='inhibitor', type_act='inhibitor',
                      plog_value=7.0, operator='=', units='nM', bioactivity_type='Ki'),
            _base_row(cmp='CHEMBL1', act='inhibitor', type_act='inhibitor',
                      plog_value=7.3, operator='=', units='nM', bioactivity_type='Ki'),
        ]
        ds = self._reg_ds(rows)
        ds._filt_diff_results(delta_diff=0.5)
        assert all(ds._pdresult['diff'] == 'norm')

    def test_non_equal_operator_marked_del(self):
        rows = [
            _base_row(cmp='CHEMBL1', act='inhibitor', type_act='inhibitor',
                      plog_value=7.0, operator='>=', units='nM', bioactivity_type='Ki'),
            _base_row(cmp='CHEMBL1', act='inhibitor', type_act='inhibitor',
                      plog_value=7.1, operator='=', units='nM', bioactivity_type='Ki'),
        ]
        ds = self._reg_ds(rows)
        ds._filt_diff_results(delta_diff=0.5)
        assert all(ds._pdresult['diff'] == 'del')


# ---------------------------------------------------------------------------
# ReggDataSet.mean_duplicate
# ---------------------------------------------------------------------------

class TestMeanDuplicate:
    def test_mean_computed_correctly(self):
        rows = [
            _base_row(cmp='CHEMBL1', act='inhibitor', type_act='inhibitor',
                      plog_value=7.0, value=100.0, units='nM', bioactivity_type='Ki'),
            _base_row(cmp='CHEMBL1', act='inhibitor', type_act='inhibitor',
                      plog_value=8.0, value=10.0, units='nM', bioactivity_type='Ki'),
        ]
        ds = _make_reg_ds(rows)
        ds.mean_duplicate()
        assert ds._pdresult['mean_plog_value'].iloc[0] == pytest.approx(7.5)
        assert ds._pdresult['mean_value'].iloc[0] == pytest.approx(55.0)

    def test_single_row_unchanged(self):
        rows = [_base_row(cmp='CHEMBL1', act='inhibitor', type_act='inhibitor',
                          plog_value=7.0, value=100.0, units='nM', bioactivity_type='Ki')]
        ds = _make_reg_ds(rows)
        ds.mean_duplicate()
        assert ds._pdresult['mean_plog_value'].iloc[0] == pytest.approx(7.0)
