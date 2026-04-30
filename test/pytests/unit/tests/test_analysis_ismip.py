"""
Unit tests for pyissm.analysis.ismip - ISMIP analysis functions.
"""

import pytest
import numpy as np
from types import SimpleNamespace

try:
    from pyissm.analysis.ismip import calc_perc_ice_cover, get_ismip_variable
    ISMIP_AVAILABLE = True
except ImportError:
    ISMIP_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not ISMIP_AVAILABLE,
    reason="pyissm.analysis.ismip not available"
)


# ============== CALC_PERC_ICE_COVER TESTS ==============

class TestCalcPercIceCover:
    """Tests for calc_perc_ice_cover()."""

    def test_full_coverage(self):
        # When ice_area == total_ice_area, should return 100%
        total = np.array([1000.0, 2000.0])
        ice = np.array([1000.0, 2000.0])
        result = calc_perc_ice_cover(total, ice)
        assert np.allclose(result, 100.0)

    def test_half_coverage(self):
        total = np.array([1000.0])
        ice = np.array([500.0])
        result = calc_perc_ice_cover(total, ice)
        assert np.isclose(result[0], 50.0)

    def test_zero_coverage(self):
        total = np.array([1000.0])
        ice = np.array([0.0])
        result = calc_perc_ice_cover(total, ice)
        assert np.isclose(result[0], 0.0)

    def test_scalar_inputs(self):
        result = calc_perc_ice_cover(1000.0, 250.0)
        assert np.isclose(result, 25.0)

    def test_array_output_shape(self):
        total = np.ones(5) * 1000.0
        ice = np.arange(5) * 200.0
        result = calc_perc_ice_cover(total, ice)
        assert result.shape == (5,)

    def test_values_are_percentage(self):
        total = np.array([100.0, 200.0, 400.0])
        ice = np.array([25.0, 50.0, 100.0])
        result = calc_perc_ice_cover(total, ice)
        assert np.allclose(result, 25.0)

    def test_returns_ndarray(self):
        result = calc_perc_ice_cover(np.array([500.0]), np.array([250.0]))
        assert isinstance(result, np.ndarray)


# ============== GET_ISMIP_VARIABLE TESTS ==============

def _make_md_with_areas(grounded, floating):
    """Helper to create a mock model with area data."""
    md = SimpleNamespace()
    md.results = SimpleNamespace()
    md.results.TransientSolution = SimpleNamespace()
    md.results.TransientSolution.GroundedArea = np.array(grounded)
    md.results.TransientSolution.FloatingArea = np.array(floating)
    return md


class TestGetIsmipVariable:
    """Tests for get_ismip_variable()."""

    def test_land_ice_area_fraction_returns_100(self):
        md = _make_md_with_areas([800.0, 900.0], [200.0, 100.0])
        result = get_ismip_variable(md, 'land_ice_area_fraction')
        assert result is not None
        assert np.allclose(result, 100.0)

    def test_floating_ice_shelf_area_fraction(self):
        md = _make_md_with_areas([800.0], [200.0])
        result = get_ismip_variable(md, 'floating_ice_shelf_area_fraction')
        assert result is not None
        assert np.isclose(result[0], 20.0)  # 200/(800+200)*100

    def test_grounded_ice_sheet_area_fraction(self):
        md = _make_md_with_areas([800.0], [200.0])
        result = get_ismip_variable(md, 'grounded_ice_sheet_area_fraction')
        assert result is not None
        assert np.isclose(result[0], 80.0)  # 800/(800+200)*100

    def test_unknown_variable_returns_none(self):
        md = _make_md_with_areas([800.0], [200.0])
        result = get_ismip_variable(md, 'unknown_variable_xyz')
        assert result is None

    def test_missing_attr_land_ice_returns_none(self, capsys):
        # Model without required attributes
        md = SimpleNamespace()
        md.results = SimpleNamespace()
        result = get_ismip_variable(md, 'land_ice_area_fraction')
        assert result is None

    def test_missing_attr_floating_returns_none(self, capsys):
        md = SimpleNamespace()
        md.results = SimpleNamespace()
        result = get_ismip_variable(md, 'floating_ice_shelf_area_fraction')
        assert result is None

    def test_missing_attr_grounded_returns_none(self, capsys):
        md = SimpleNamespace()
        md.results = SimpleNamespace()
        result = get_ismip_variable(md, 'grounded_ice_sheet_area_fraction')
        assert result is None

    def test_floating_plus_grounded_equals_100(self):
        md = _make_md_with_areas([750.0], [250.0])
        floating = get_ismip_variable(md, 'floating_ice_shelf_area_fraction')
        grounded = get_ismip_variable(md, 'grounded_ice_sheet_area_fraction')
        assert np.isclose(floating[0] + grounded[0], 100.0)

    def test_land_ice_returns_same_as_total(self):
        grounded = np.array([600.0, 700.0])
        floating = np.array([400.0, 300.0])
        md = _make_md_with_areas(grounded.tolist(), floating.tolist())
        result = get_ismip_variable(md, 'land_ice_area_fraction')
        assert np.allclose(result, 100.0)

    def test_array_inputs(self):
        md = _make_md_with_areas([500.0, 600.0, 700.0], [500.0, 400.0, 300.0])
        result = get_ismip_variable(md, 'grounded_ice_sheet_area_fraction')
        assert result is not None
        expected = np.array([50.0, 60.0, 70.0])
        assert np.allclose(result, expected)
