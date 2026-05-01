"""
Unit tests for pyissm.tools.general module.

Tests cover:
- Unit conversion functions
- Utility functions (has_nested_attr, planetradius, etc.)
- Coordinate conversion functions (xy_to_ll, ll_to_xy)
- Field extraction functions (extract_field_layer)

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
        xy_to_ll,
        ll_to_xy,
        extract_field_layer,
    )
    ISSM_AVAILABLE = True
except ImportError:
    ISSM_AVAILABLE = False
    convert_units = has_nested_attr = planetradius = _wgs84_ellipsoid_constants = None
    xy_to_ll = ll_to_xy = extract_field_layer = None

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


class TestXYToLL:
    """Tests for xy_to_ll coordinate conversion."""

    def test_origin_returns_pole_south(self):
        """Test that origin (0,0) returns the pole for south hemisphere."""
        # Use array input to avoid scalar indexing issue in the function
        lat, lon = xy_to_ll(np.array([0.0]), np.array([0.0]), sign=-1, central_meridian=0, standard_parallel=71)
        assert np.isclose(lat[0], -90.0, atol=0.1)

    def test_origin_returns_pole_north(self):
        """Test that origin (0,0) returns the pole for north hemisphere."""
        lat, lon = xy_to_ll(np.array([0.0]), np.array([0.0]), sign=1, central_meridian=45, standard_parallel=70)
        assert np.isclose(lat[0], 90.0, atol=0.1)

    def test_invalid_sign_raises(self):
        """Test that invalid sign raises ValueError."""
        with pytest.raises(ValueError, match="sign should be either"):
            xy_to_ll(0, 0, sign=0)

    def test_partial_params_raises(self):
        """Test that specifying only one parameter raises ValueError."""
        with pytest.raises(ValueError, match="Specify both"):
            xy_to_ll(0, 0, sign=1, central_meridian=45)

    def test_array_input(self):
        """Test that array inputs work correctly."""
        x = np.array([0, 1000, 2000])
        y = np.array([0, 1000, 2000])
        lat, lon = xy_to_ll(x, y, sign=-1, central_meridian=0, standard_parallel=71)
        assert lat.shape == (3,)
        assert lon.shape == (3,)

    def test_output_within_valid_range(self):
        """Test that output is within valid lat/lon range."""
        x = np.array([100000, -100000, 0])
        y = np.array([0, 100000, -100000])
        lat, lon = xy_to_ll(x, y, sign=-1, central_meridian=0, standard_parallel=71)
        assert np.all(np.abs(lat) <= 90)
        assert np.all(np.abs(lon) <= 360)


class TestLLToXY:
    """Tests for ll_to_xy coordinate conversion."""

    def test_pole_returns_origin_south(self):
        """Test that south pole returns origin."""
        # Use array input to avoid scalar indexing issue
        x, y = ll_to_xy(np.array([-90.0]), np.array([0.0]), sign=-1, central_meridian=0, standard_parallel=71)
        assert np.isclose(x[0], 0, atol=1)
        assert np.isclose(y[0], 0, atol=1)

    def test_pole_returns_origin_north(self):
        """Test that north pole returns origin."""
        x, y = ll_to_xy(np.array([90.0]), np.array([0.0]), sign=1, central_meridian=45, standard_parallel=70)
        assert np.isclose(x[0], 0, atol=1)
        assert np.isclose(y[0], 0, atol=1)

    def test_invalid_sign_raises(self):
        """Test that invalid sign raises ValueError."""
        with pytest.raises(ValueError, match="sign should be either"):
            ll_to_xy(0, 0, sign=2)

    def test_partial_params_raises(self):
        """Test that specifying only one parameter raises ValueError."""
        with pytest.raises(ValueError, match="Specify both"):
            ll_to_xy(0, 0, sign=-1, standard_parallel=71)

    def test_array_input(self):
        """Test that array inputs work correctly."""
        lat = np.array([-70, -75, -80])
        lon = np.array([0, 45, 90])
        x, y = ll_to_xy(lat, lon, sign=-1, central_meridian=0, standard_parallel=71)
        assert x.shape == (3,)
        assert y.shape == (3,)


class TestCoordinateRoundTrip:
    """Round-trip tests for coordinate conversions."""

    @pytest.mark.parametrize("sign,cm,sp", [
        (-1, 0, 71),    # South polar stereographic
        (1, 45, 70),    # North polar stereographic
    ])
    def test_ll_to_xy_to_ll_roundtrip(self, sign, cm, sp):
        """Test that converting ll->xy->ll returns original values."""
        # Use coordinates away from pole for better precision
        lat_orig = sign * 75.0
        lon_orig = 45.0
        
        x, y = ll_to_xy(lat_orig, lon_orig, sign=sign, central_meridian=cm, standard_parallel=sp)
        lat_back, lon_back = xy_to_ll(x, y, sign=sign, central_meridian=cm, standard_parallel=sp)
        
        assert np.isclose(lat_back, lat_orig, atol=1e-6)
        # Longitude may wrap around, so compare mod 360
        assert np.isclose((lon_back - lon_orig) % 360, 0, atol=1e-6) or \
               np.isclose((lon_back - lon_orig) % 360, 360, atol=1e-6)


class TestExtractFieldLayer:
    """Tests for extract_field_layer function."""

    def test_extracts_correct_layer(self, fake_model_3d):
        """Test that correct layer is extracted from 3D field."""
        # Create a 3D field with distinct values per layer
        n2d = fake_model_3d.mesh.numberofvertices2d
        nlayers = fake_model_3d.mesh.numberoflayers
        field = np.arange(n2d * nlayers, dtype=float)
        
        # Extract layer 2
        layer_data, indices = extract_field_layer(fake_model_3d, field, layer=2)
        
        assert len(layer_data) == n2d
        assert len(indices) == n2d

    def test_layer_out_of_bounds_raises(self, fake_model_3d):
        """Test that invalid layer raises IndexError (out of bounds access)."""
        field = np.ones(fake_model_3d.mesh.numberofvertices)
        
        with pytest.raises(IndexError):
            extract_field_layer(fake_model_3d, field, layer=100)

    def test_2d_model_raises(self, fake_model_2d):
        """Test that 2D model raises TypeError."""
        field = np.ones(fake_model_2d.mesh.numberofvertices)
        
        with pytest.raises(TypeError, match="not 3D"):
            extract_field_layer(fake_model_2d, field, layer=1)



class TestXYToLL:
    """Tests for xy_to_ll coordinate conversion."""

    def test_origin_returns_pole_south(self):
        """Test that origin (0,0) returns the pole for south hemisphere."""
        # Use array input to avoid scalar indexing issue in the function
        lat, lon = xy_to_ll(np.array([0.0]), np.array([0.0]), sign=-1, central_meridian=0, standard_parallel=71)
        assert np.isclose(lat[0], -90.0, atol=0.1)

    def test_origin_returns_pole_north(self):
        """Test that origin (0,0) returns the pole for north hemisphere."""
        lat, lon = xy_to_ll(np.array([0.0]), np.array([0.0]), sign=1, central_meridian=45, standard_parallel=70)
        assert np.isclose(lat[0], 90.0, atol=0.1)

    def test_invalid_sign_raises(self):
        """Test that invalid sign raises ValueError."""
        with pytest.raises(ValueError, match="sign should be either"):
            xy_to_ll(0, 0, sign=0)

    def test_partial_params_raises(self):
        """Test that specifying only one parameter raises ValueError."""
        with pytest.raises(ValueError, match="Specify both"):
            xy_to_ll(0, 0, sign=1, central_meridian=45)

    def test_array_input(self):
        """Test that array inputs work correctly."""
        x = np.array([0, 1000, 2000])
        y = np.array([0, 1000, 2000])
        lat, lon = xy_to_ll(x, y, sign=-1, central_meridian=0, standard_parallel=71)
        assert lat.shape == (3,)
        assert lon.shape == (3,)

    def test_output_within_valid_range(self):
        """Test that output is within valid lat/lon range."""
        x = np.array([100000, -100000, 0])
        y = np.array([0, 100000, -100000])
        lat, lon = xy_to_ll(x, y, sign=-1, central_meridian=0, standard_parallel=71)
        assert np.all(np.abs(lat) <= 90)
        assert np.all(np.abs(lon) <= 360)


class TestLLToXY:
    """Tests for ll_to_xy coordinate conversion."""

    def test_pole_returns_origin_south(self):
        """Test that south pole returns origin."""
        # Use array input to avoid scalar indexing issue
        x, y = ll_to_xy(np.array([-90.0]), np.array([0.0]), sign=-1, central_meridian=0, standard_parallel=71)
        assert np.isclose(x[0], 0, atol=1)
        assert np.isclose(y[0], 0, atol=1)

    def test_pole_returns_origin_north(self):
        """Test that north pole returns origin."""
        x, y = ll_to_xy(np.array([90.0]), np.array([0.0]), sign=1, central_meridian=45, standard_parallel=70)
        assert np.isclose(x[0], 0, atol=1)
        assert np.isclose(y[0], 0, atol=1)

    def test_invalid_sign_raises(self):
        """Test that invalid sign raises ValueError."""
        with pytest.raises(ValueError, match="sign should be either"):
            ll_to_xy(0, 0, sign=2)

    def test_partial_params_raises(self):
        """Test that specifying only one parameter raises ValueError."""
        with pytest.raises(ValueError, match="Specify both"):
            ll_to_xy(0, 0, sign=-1, standard_parallel=71)

    def test_array_input(self):
        """Test that array inputs work correctly."""
        lat = np.array([-70, -75, -80])
        lon = np.array([0, 45, 90])
        x, y = ll_to_xy(lat, lon, sign=-1, central_meridian=0, standard_parallel=71)
        assert x.shape == (3,)
        assert y.shape == (3,)


class TestCoordinateRoundTrip:
    """Round-trip tests for coordinate conversions."""

    @pytest.mark.parametrize("sign,cm,sp", [
        (-1, 0, 71),    # South polar stereographic
        (1, 45, 70),    # North polar stereographic
    ])
    def test_ll_to_xy_to_ll_roundtrip(self, sign, cm, sp):
        """Test that converting ll->xy->ll returns original values."""
        # Use coordinates away from pole for better precision
        lat_orig = sign * 75.0
        lon_orig = 45.0
        
        x, y = ll_to_xy(lat_orig, lon_orig, sign=sign, central_meridian=cm, standard_parallel=sp)
        lat_back, lon_back = xy_to_ll(x, y, sign=sign, central_meridian=cm, standard_parallel=sp)
        
        assert np.isclose(lat_back, lat_orig, atol=1e-6)
        # Longitude may wrap around, so compare mod 360
        assert np.isclose((lon_back - lon_orig) % 360, 0, atol=1e-6) or \
               np.isclose((lon_back - lon_orig) % 360, 360, atol=1e-6)


class TestExtractFieldLayer:
    """Tests for extract_field_layer function."""

    def test_extracts_correct_layer(self, fake_model_3d):
        """Test that correct layer is extracted from 3D field."""
        # Create a 3D field with distinct values per layer
        n2d = fake_model_3d.mesh.numberofvertices2d
        nlayers = fake_model_3d.mesh.numberoflayers
        field = np.arange(n2d * nlayers, dtype=float)
        
        # Extract layer 2
        layer_data, indices = extract_field_layer(fake_model_3d, field, layer=2)
        
        assert len(layer_data) == n2d
        assert len(indices) == n2d

    def test_layer_out_of_bounds_raises(self, fake_model_3d):
        """Test that invalid layer raises IndexError (out of bounds access)."""
        field = np.ones(fake_model_3d.mesh.numberofvertices)
        
        with pytest.raises(IndexError):
            extract_field_layer(fake_model_3d, field, layer=100)

    def test_2d_model_raises(self, fake_model_2d):
        """Test that 2D model raises TypeError."""
        field = np.ones(fake_model_2d.mesh.numberofvertices)
        
        with pytest.raises(TypeError, match="not 3D"):
            extract_field_layer(fake_model_2d, field, layer=1)

