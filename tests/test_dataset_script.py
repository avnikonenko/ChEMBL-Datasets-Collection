import pytest
from dataset_script import cuttof_target_type


class TestCutoffTargetType:
    def test_epigenetic_factor(self):
        act, inact, type_act, contr = cuttof_target_type('epigenetic_factor')
        assert act['Ki'] == (6, '>=')
        assert 'inhibitor' in type_act

    def test_transporter(self):
        act, inact, type_act, contr = cuttof_target_type('transporter')
        assert act['Ki'] == (6, '>=')

    def test_ion_channel(self):
        act, inact, type_act, contr = cuttof_target_type('ion_channel')
        assert act['Ki'] == (5, '>=')
        assert 'blockator' in type_act

    def test_enzyme(self):
        act, inact, type_act, contr = cuttof_target_type('enzyme')
        assert act['Ki'] == (7.5, '>=')
        assert type_act == ['allosteric', 'activator', 'inhibitor', 'binding']
        assert contr_list_has(contr, 'activator', 'inhibitor')

    def test_membrane_protein(self):
        act, inact, type_act, contr = cuttof_target_type('membrane_protein')
        assert act['Ki'] == (7, '>=')

    def test_membrane_receptor(self):
        act, inact, type_act, contr = cuttof_target_type('membrane_receptor')
        assert act['Ki'] == (7, '>=')

    def test_unknown_class_uses_default(self):
        # previously dead code due to the 'membrane_protein' or ... bug
        act, inact, type_act, contr = cuttof_target_type('nuclear_receptor')
        assert act['Ki'] == (6, '>=')

    def test_space_normalization(self):
        act1, _, _, _ = cuttof_target_type('ion channel')
        act2, _, _, _ = cuttof_target_type('ion_channel')
        assert act1 == act2

    def test_case_normalization(self):
        act1, _, _, _ = cuttof_target_type('Enzyme')
        act2, _, _, _ = cuttof_target_type('enzyme')
        assert act1 == act2

    def test_enzyme_not_matched_by_membrane_branch(self):
        # regression: before the fix, enzyme would fall into the membrane branch
        # because 'membrane_protein' or ... was always True
        _, _, type_act, _ = cuttof_target_type('enzyme')
        assert 'blockator' not in type_act  # blockator belongs to ion_channel
        assert 'activator' in type_act      # activator belongs to enzyme


def contr_list_has(contr_list, a, b):
    return any((x == a and y == b) or (x == b and y == a) for x, y in contr_list)
