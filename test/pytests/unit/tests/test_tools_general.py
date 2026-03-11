"""
Unit tests for pyissm.tools.general module.

Tests cover:
- Unit conversion functions
- Utility functions (has_nested_attr, planetradius, etc.)

Note: These tests require the ISSM backend to be available.
"""

import numpy as np
import pytest
from types import SimpleNamespace

try:
    from pyissm.tools.general import (
        convert_units,
        has_nested_attr,
        planetradius,
        _wgs84_ellipsoid_constants,
    )
    ISSM_AVAILABLE = True
except ImportError:
    ISSM_AVAILABLE = False
    convert_units = has_nested_attr = planetradius = _wgs84_ellipsoid_constants = None

pytestmark = pytest.mark.skipif(
    not ISSM_AVAILABLE,
    reason="ISSM backend not available"
)


class TestConvertUnits:
    """Tests for the convert_units function."""

    # ---------------------
    # Length conversions
    # ---------------------
    def test_m_to_km(self):
        """Test meters to kilometers conversion."""
        result = convert_units('m', 'km', 1000.0)
        assert np.isclose(result, 1.0)

    def test_km_to_m(self):
        """Test kilometers to meters conversion."""
        result = convert_units('km', 'm', 1.0)
        assert np.isclose(result, 1000.0)

    def test_m_to_km_array(self):
        """Test meters to kilometers with array input."""
        data = np.array([1000.0, 2000.0, 5000.0])
        result = convert_units('m', 'km', data)
        expected = np.array([1.0, 2.0, 5.0])
        np.testing.assert_array_almost_equal(result, expected)

    # ---------------------
    # Area conversions
    # ---------------------
    def test_m2_to_km2(self):
        """Test square meters to square kilometers conversion."""
        result = convert_units('m2', 'km2', 1e6)
        assert np.isclose(result, 1.0)

    def test_km2_to_m2(self):
        """Test square kilometers to square meters conversion."""
        result = convert_units('km2', 'm2', 1.0)
        assert np.isclose(result, 1e6)

    # ---------------------
    # Speed conversions
    # ---------------------
    def test_ms_to_myr(self):
        """Test m/s to m/yr conversion."""
        yts = 365 * 24 * 60 * 60  # seconds per year
        result = convert_units('ms-1', 'myr-1', 1.0)
        assert np.isclose(result, yts)

    def test_myr_to_ms(self):
        """Test m/yr to m/s conversion."""
        yts = 365 * 24 * 60 * 60
        result = convert_units('myr-1', 'ms-1', yts)
        assert np.isclose(result, 1.0)

    # ---------------------
    # Volume conversions
    # ---------------------
    def test_m3_to_km3(self):
        """Test cubic meters to cubic kilometers conversion."""
        result = convert_units('m3', 'km3', 1e9)
        assert np.isclose(result, 1.0)

    def test_km3_to_m3(self):
        """Test cubic kilometers to cubic meters conversion."""
        result = convert_units('km3', 'm3', 1.0)
        assert np.isclose(result, 1e9)

    # ---------------------
    # Mass conversions
    # ---------------------
    def test_gt_to_km3(self):
        """Test Gigatons to cubic kilometers conversion (ice density)."""
        # 1 Gt = 1e12 kg, with rho_ice = 917 kg/m3
        # 1 Gt -> 1e12 / 917 m3 -> (1e12 / 917) / 1e9 km3
        result = convert_units('Gt', 'km3', 1.0)
        expected = (1e12 / 917) / 1e9
        assert np.isclose(result, expected)

    def test_km3_to_gt(self):
        """Test cubic kilometers to Gigatons conversion."""
        # Inverse of above
        gt_value = 1.0
        km3_value = convert_units('Gt', 'km3', gt_value)
        result = convert_units('km3', 'Gt', km3_value)
        assert np.isclose(result, gt_value)

    def test_m3_to_kg(self):
        """Test cubic meters to kilograms (ice)."""
        result = convert_units('m3', 'kg', 1.0, rho_ice=917)
        assert np.isclose(result, 917)

    def test_kg_to_m3(self):
        """Test kilograms to cubic meters (ice)."""
        result = convert_units('kg', 'm3', 917, rho_ice=917)
        assert np.isclose(result, 1.0)

    # ---------------------
    # Rate conversions
    # ---------------------
    def test_gtyr_to_kgs(self):
        """Test Gt/yr to kg/s conversion."""
        yts = 365 * 24 * 60 * 60
        result = convert_units('Gtyr-1', 'kgs-1', 1.0)
        expected = 1e12 / yts
        assert np.isclose(result, expected)

    def test_kgs_to_gtyr(self):
        """Test kg/s to Gt/yr conversion."""
        yts = 365 * 24 * 60 * 60
        kgs = 1e12 / yts  # 1 Gt/yr in kg/s
        result = convert_units('kgs-1', 'Gtyr-1', kgs)
        assert np.isclose(result, 1.0)

    # ---------------------
    # Error handling
    # ---------------------
    def test_invalid_input_units_raises(self):
        """Test that invalid input units raise ValueError."""
        with pytest.raises(ValueError, match="Invalid input_units"):
            convert_units('invalid', 'km', 1.0)

    def test_invalid_output_units_raises(self):
        """Test that invalid output units raise ValueError."""
        with pytest.raises(ValueError, match="Invalid output_units"):
            convert_units('m', 'invalid', 1.0)

    def test_unsupported_conversion_raises(self):
        """Test that unsupported conversion raises ValueError."""
        with pytest.raises(ValueError, match="not supported"):
            convert_units('m', 'Gt', 1.0)

    # ---------------------
    # Custom parameters
    # ---------------------
    def test_custom_yts(self):
        """Test conversion with custom seconds-per-year."""
        custom_yts = 365.25 * 24 * 60 * 60  # Include leap year
        result = convert_units('ms-1', 'myr-1', 1.0, yts=custom_yts)
        assert np.isclose(result, custom_yts)

    def test_custom_rho_ice(self):
        """Test conversion with custom ice density."""
        custom_rho = 910  # Different ice density
        result = convert_units('m3', 'kg', 1.0, rho_ice=custom_rho)
        assert np.isclose(result, custom_rho)


class TestHasNestedAttr:
    """Tests for the has_nested_attr function."""

    def test_single_attribute_exists(self):
        """Test checking a single existing attribute."""
        obj = SimpleNamespace(a=1)
        assert has_nested_attr(obj, 'a') is True

    def test_single_attribute_missing(self):
        """Test checking a single missing attribute."""
        obj = SimpleNamespace(a=1)
        assert has_nested_attr(obj, 'b') is False

    def test_nested_attribute_exists(self):
        """Test checking nested existing attributes."""
        obj = SimpleNamespace(a=SimpleNamespace(b=SimpleNamespace(c=1)))
        assert has_nested_attr(obj, 'a', 'b', 'c') is True

    def test_nested_attribute_partial_exists(self):
        """Test when only part of the chain exists."""
        obj = SimpleNamespace(a=SimpleNamespace(b=1))
        assert has_nested_attr(obj, 'a', 'b', 'c') is False

    def test_nested_attribute_missing_middle(self):
        """Test when middle attribute is missing."""
        obj = SimpleNamespace(a=1)
        assert has_nested_attr(obj, 'a', 'b', 'c') is False

    def test_empty_attrs(self):
        """Test with no attributes to check."""
        obj = SimpleNamespace(a=1)
        assert has_nested_attr(obj) is True

    def test_with_model_like_structure(self, fake_model_2d):
        """Test with a model-like nested structure."""
        assert has_nested_attr(fake_model_2d, 'mesh', 'x') is True
        assert has_nested_attr(fake_model_2d, 'mesh', 'elements2d') is False
        assert has_nested_attr(fake_model_2d, 'geometry', 'surface') is True


class TestPlanetRadius:
    """Tests for the planetradius function."""

    def test_earth_radius(self):
        """Test Earth's radius."""
        radius = planetradius('earth')
        assert np.isclose(radius, 6.371012e6)

    def test_europa_radius(self):
        """Test Europa's radius."""
        radius = planetradius('europa')
        assert np.isclose(radius, 1.5008e6)

    def test_invalid_planet_raises(self):
        """Test that invalid planet raises TypeError."""
        with pytest.raises(TypeError, match="not supported"):
            planetradius('mars')


class TestWGS84Constants:
    """Tests for WGS84 ellipsoid constants."""

    def test_equatorial_radius(self):
        """Test WGS84 equatorial radius."""
        re, f, ex2, ex = _wgs84_ellipsoid_constants()
        assert re == 6378137

    def test_flattening(self):
        """Test WGS84 flattening."""
        re, f, ex2, ex = _wgs84_ellipsoid_constants()
        expected_f = 1.0 / 298.257223563
        assert np.isclose(f, expected_f)

    def test_eccentricity_squared(self):
        """Test WGS84 eccentricity squared."""
        re, f, ex2, ex = _wgs84_ellipsoid_constants()
        expected_ex2 = 2 * f - f**2
        assert np.isclose(ex2, expected_ex2)

    def test_eccentricity(self):
        """Test WGS84 eccentricity."""
        re, f, ex2, ex = _wgs84_ellipsoid_constants()
        assert np.isclose(ex, np.sqrt(ex2))


class TestConvertUnitsRoundTrip:
    """Round-trip tests to verify conversion consistency."""

    @pytest.mark.parametrize("units_pair", [
        ('m', 'km'),
        ('m2', 'km2'),
        ('ms-1', 'myr-1'),
        ('m3', 'km3'),
        ('Gt', 'km3'),
        ('Gtyr-1', 'kgs-1'),
        ('m3', 'kg'),
    ])
    def test_round_trip_conversion(self, units_pair):
        """Test that converting and back returns original value."""
        unit_a, unit_b = units_pair
        original = 42.0
        
        converted = convert_units(unit_a, unit_b, original)
        recovered = convert_units(unit_b, unit_a, converted)
        
        assert np.isclose(recovered, original, rtol=1e-10)
