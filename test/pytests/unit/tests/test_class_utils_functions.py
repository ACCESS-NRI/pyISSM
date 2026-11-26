"""
Unit tests for pyissm.model.classes.class_utils - pure utility functions.
"""

import pytest
import numpy as np
from types import SimpleNamespace

try:
    from pyissm.model.classes.class_utils import (
        marshall_inversion_cost_functions,
        supported_inversion_control_parameters,
        supported_inversion_cost_functions,
        supported_analyses,
        supported_stochastic_forcings,
        _check_values,
        _check_bound,
        _check_size,
    )
    CLASS_UTILS_AVAILABLE = True
except ImportError:
    CLASS_UTILS_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not CLASS_UTILS_AVAILABLE,
    reason="pyissm.model.classes.class_utils not available"
)


def _make_mock_md(nv=10, ne=8):
    """Helper: create a minimal mock model for check functions."""
    md = SimpleNamespace()
    md.mesh = SimpleNamespace()
    md.mesh.numberofvertices = nv
    md.mesh.numberofelements = ne
    md.private = SimpleNamespace()
    md.private.isconsistent = True
    messages = []
    def check_message(s):
        md.private.isconsistent = False
        messages.append(s)
    md.check_message = check_message
    md._messages = messages
    return md


# ============== MARSHALL_INVERSION_COST_FUNCTIONS ==============

class TestMarshallInversionCostFunctions:
    """Tests for marshall_inversion_cost_functions()."""

    def test_single_int_101(self):
        result = marshall_inversion_cost_functions(101)
        assert result == ['SurfaceAbsVelMisfit']

    def test_single_int_201(self):
        result = marshall_inversion_cost_functions(201)
        assert result == ['ThicknessAbsMisfit']

    def test_single_int_501(self):
        result = marshall_inversion_cost_functions(501)
        assert result == ['DragCoefficientAbsGradient']

    def test_single_int_502(self):
        result = marshall_inversion_cost_functions(502)
        assert result == ['RheologyBbarAbsGradient']

    def test_single_int_503(self):
        result = marshall_inversion_cost_functions(503)
        assert result == ['ThicknessAbsGradient']

    def test_list_of_codes(self):
        result = marshall_inversion_cost_functions([101, 201, 501])
        assert result == ['SurfaceAbsVelMisfit', 'ThicknessAbsMisfit', 'DragCoefficientAbsGradient']

    def test_returns_list(self):
        result = marshall_inversion_cost_functions(101)
        assert isinstance(result, list)

    def test_102_surfacerelvel(self):
        result = marshall_inversion_cost_functions(102)
        assert result == ['SurfaceRelVelMisfit']

    def test_103_surfacelogvel(self):
        result = marshall_inversion_cost_functions(103)
        assert result == ['SurfaceLogVelMisfit']

    def test_104_surfacelogvxvy(self):
        result = marshall_inversion_cost_functions(104)
        assert result == ['SurfaceLogVxVyMisfit']

    def test_105_surfaceaveravel(self):
        result = marshall_inversion_cost_functions(105)
        assert result == ['SurfaceAverageVelMisfit']

    def test_504_505(self):
        result = marshall_inversion_cost_functions([504, 505])
        assert result == ['ThicknessAlongGradient', 'ThicknessAcrossGradient']

    def test_empty_list(self):
        result = marshall_inversion_cost_functions([])
        assert result == []


# ============== SUPPORTED_INVERSION_CONTROL_PARAMETERS ==============

class TestSupportedInversionControlParameters:
    """Tests for supported_inversion_control_parameters()."""

    def test_returns_list(self):
        result = supported_inversion_control_parameters()
        assert isinstance(result, list)

    def test_non_empty(self):
        result = supported_inversion_control_parameters()
        assert len(result) > 0

    def test_all_strings(self):
        result = supported_inversion_control_parameters()
        assert all(isinstance(s, str) for s in result)

    def test_contains_friction_coefficient(self):
        result = supported_inversion_control_parameters()
        assert 'FrictionCoefficient' in result

    def test_contains_materials_rheology_bbar(self):
        result = supported_inversion_control_parameters()
        assert 'MaterialsRheologyBbar' in result

    def test_contains_thickness(self):
        result = supported_inversion_control_parameters()
        assert 'Thickness' in result


# ============== SUPPORTED_INVERSION_COST_FUNCTIONS ==============

class TestSupportedInversionCostFunctions:
    """Tests for supported_inversion_cost_functions()."""

    def test_returns_list(self):
        result = supported_inversion_cost_functions()
        assert isinstance(result, list)

    def test_non_empty(self):
        result = supported_inversion_cost_functions()
        assert len(result) > 0

    def test_all_ints(self):
        result = supported_inversion_cost_functions()
        assert all(isinstance(c, int) for c in result)

    def test_contains_101(self):
        result = supported_inversion_cost_functions()
        assert 101 in result

    def test_contains_201(self):
        result = supported_inversion_cost_functions()
        assert 201 in result

    def test_contains_501_to_508(self):
        result = supported_inversion_cost_functions()
        for code in range(501, 509):
            assert code in result

    def test_contains_601_not_602_to_604(self):
        # 602-604 are in upstream's supportedcostfunctions but have no name to marshall (issue #78)
        result = supported_inversion_cost_functions()
        assert 601 in result
        for code in range(602, 605):
            assert code not in result


# ============== SUPPORTED_ANALYSES ==============

class TestSupportedAnalyses:
    """Tests for supported_analyses()."""

    def test_returns_list(self):
        result = supported_analyses()
        assert isinstance(result, list)

    def test_non_empty(self):
        result = supported_analyses()
        assert len(result) > 0

    def test_all_strings(self):
        result = supported_analyses()
        assert all(isinstance(s, str) for s in result)

    def test_contains_stressbalance(self):
        result = supported_analyses()
        assert 'StressbalanceAnalysis' in result

    def test_contains_masstransport(self):
        result = supported_analyses()
        assert 'MasstransportAnalysis' in result

    def test_contains_thermal(self):
        result = supported_analyses()
        assert 'ThermalAnalysis' in result

    def test_contains_levelset(self):
        result = supported_analyses()
        assert 'LevelsetAnalysis' in result


# ============== SUPPORTED_STOCHASTIC_FORCINGS ==============

class TestSupportedStochasticForcings:
    """Tests for supported_stochastic_forcings()."""

    def test_default_returns_list(self):
        result = supported_stochastic_forcings()
        assert isinstance(result, list)

    def test_return_dict_true_returns_dict(self):
        result = supported_stochastic_forcings(return_dict=True)
        assert isinstance(result, dict)

    def test_list_contains_strings(self):
        result = supported_stochastic_forcings()
        assert all(isinstance(s, str) for s in result)

    def test_contains_smbarma(self):
        result = supported_stochastic_forcings()
        assert 'SMBarma' in result

    def test_contains_smbforcing(self):
        result = supported_stochastic_forcings()
        assert 'SMBforcing' in result

    def test_dict_keys_match_list(self):
        keys = supported_stochastic_forcings(return_dict=True).keys()
        lst = supported_stochastic_forcings()
        assert set(keys) == set(lst)


# ============== _CHECK_VALUES ==============

class TestCheckValues:
    """Tests for _check_values()."""

    def test_valid_value_passes(self):
        md = _make_mock_md()
        _check_values(md, 0, 'test.field', [0, 1])
        assert md.private.isconsistent is True

    def test_invalid_value_fails(self):
        md = _make_mock_md()
        _check_values(md, 5, 'test.field', [0, 1])
        assert md.private.isconsistent is False

    def test_array_all_valid(self):
        md = _make_mock_md()
        _check_values(md, np.array([0, 1, 0, 1]), 'test.field', [0, 1])
        assert md.private.isconsistent is True

    def test_array_some_invalid(self):
        md = _make_mock_md()
        _check_values(md, np.array([0, 1, 2]), 'test.field', [0, 1])
        assert md.private.isconsistent is False

    def test_custom_message(self):
        md = _make_mock_md()
        _check_values(md, 99, 'test.field', [0, 1], message='Custom error')
        assert 'Custom error' in md._messages[0]


# ============== _CHECK_BOUND ==============

class TestCheckBound:
    """Tests for _check_bound()."""

    def test_ge_valid(self):
        md = _make_mock_md()
        _check_bound(md, 5.0, 'test.field', '>=', 0.0, True)
        assert md.private.isconsistent is True

    def test_ge_invalid(self):
        md = _make_mock_md()
        _check_bound(md, -1.0, 'test.field', '>=', 0.0, False)
        assert md.private.isconsistent is False

    def test_gt_valid(self):
        md = _make_mock_md()
        _check_bound(md, 1.0, 'test.field', '>', 0.0, False)
        assert md.private.isconsistent is True

    def test_gt_invalid(self):
        md = _make_mock_md()
        _check_bound(md, 0.0, 'test.field', '>', 0.0, False)
        assert md.private.isconsistent is False

    def test_le_valid(self):
        md = _make_mock_md()
        _check_bound(md, 5.0, 'test.field', '<=', 10.0, False)
        assert md.private.isconsistent is True

    def test_le_invalid(self):
        md = _make_mock_md()
        _check_bound(md, 15.0, 'test.field', '<=', 10.0, False)
        assert md.private.isconsistent is False

    def test_lt_valid(self):
        md = _make_mock_md()
        _check_bound(md, 5.0, 'test.field', '<', 10.0, False)
        assert md.private.isconsistent is True

    def test_lt_invalid(self):
        md = _make_mock_md()
        _check_bound(md, 10.0, 'test.field', '<', 10.0, False)
        assert md.private.isconsistent is False

    def test_nan_with_allow_nan_true(self):
        md = _make_mock_md()
        # NaN values should be skipped when allow_nan=True
        _check_bound(md, np.nan, 'test.field', '>=', 0.0, True)
        assert md.private.isconsistent is True

    def test_array_all_valid(self):
        md = _make_mock_md()
        _check_bound(md, np.array([1.0, 2.0, 3.0]), 'test.field', '>=', 0.0, False)
        assert md.private.isconsistent is True

    def test_array_some_invalid(self):
        md = _make_mock_md()
        _check_bound(md, np.array([1.0, -1.0, 3.0]), 'test.field', '>=', 0.0, False)
        assert md.private.isconsistent is False


# ============== _CHECK_SIZE ==============

class TestCheckSize:
    """Tests for _check_size()."""

    def test_correct_size_passes(self):
        md = _make_mock_md()
        field = np.zeros((3, 4))
        _check_size(md, field, 'test.field', (3, 4))
        assert md.private.isconsistent is True

    def test_wrong_size_fails(self):
        md = _make_mock_md()
        field = np.zeros((3, 4))
        _check_size(md, field, 'test.field', (4, 3))
        assert md.private.isconsistent is False

    def test_none_expected_passes(self):
        md = _make_mock_md()
        field = np.zeros((3, 4))
        _check_size(md, field, 'test.field', None)
        assert md.private.isconsistent is True

    def test_nan_wildcard_in_expected(self):
        md = _make_mock_md()
        field = np.zeros((3, 4))
        _check_size(md, field, 'test.field', (3, np.nan))
        assert md.private.isconsistent is True

    def test_wrong_ndim_fails(self):
        md = _make_mock_md()
        field = np.zeros((3,))
        _check_size(md, field, 'test.field', (3, 4))
        assert md.private.isconsistent is False


class TestInversionCostFunctionsConsistent:
    """Every cost function accepted by check_consistency can be marshalled (issue #78)."""

    def test_every_supported_code_marshalls(self):
        codes = supported_inversion_cost_functions()
        assert len(marshall_inversion_cost_functions(codes)) == len(codes)

    def test_dict_matches_list(self):
        mapping = supported_inversion_cost_functions(return_dict = True)
        assert list(mapping) == supported_inversion_cost_functions()

    @pytest.mark.parametrize('code, name', [(506, 'BalancethicknessMisfit'),
                                            (507, 'RheologyBAbsGradient'),
                                            (508, 'RheologyBInitialguessMisfit'),
                                            (510, 'ThicknessPositive'),
                                            (601, 'SurfaceAbsMisfit')])
    def test_codes_previously_missing(self, code, name):
        assert marshall_inversion_cost_functions(code) == [name]

    def test_numpy_input(self):
        import numpy as np
        assert marshall_inversion_cost_functions(np.array([101., 501.])) == ['SurfaceAbsVelMisfit', 'DragCoefficientAbsGradient']
        assert marshall_inversion_cost_functions(np.int64(201)) == ['ThicknessAbsMisfit']

    @pytest.mark.parametrize('code', [602, 999])
    def test_unsupported_code_raises(self, code):
        assert code not in supported_inversion_cost_functions()
        with pytest.raises(ValueError, match = 'unsupported cost function'):
            marshall_inversion_cost_functions(code)
