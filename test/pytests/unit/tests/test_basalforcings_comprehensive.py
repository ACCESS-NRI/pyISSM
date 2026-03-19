"""
Comprehensive unit tests for pyissm.model.classes.basalforcings module.

Tests cover all basalforcings classes.
"""

import pytest

try:
    from pyissm.model.classes import basalforcings
    BF_AVAILABLE = True
except ImportError:
    BF_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not BF_AVAILABLE,
    reason="Basalforcings classes not available"
)


class TestBasalforcingsDefaultComprehensive:
    """Comprehensive tests for default basalforcings class."""

    def test_init(self):
        b = basalforcings.default()
        assert b is not None

    def test_has_geothermalflux(self):
        b = basalforcings.default()
        assert hasattr(b, 'geothermalflux')

    def test_has_groundedice_melting_rate(self):
        b = basalforcings.default()
        assert hasattr(b, 'groundedice_melting_rate')

    def test_has_floatingice_melting_rate(self):
        b = basalforcings.default()
        assert hasattr(b, 'floatingice_melting_rate')

    def test_repr(self):
        b = basalforcings.default()
        r = repr(b)
        assert isinstance(r, str)

    def test_str(self):
        b = basalforcings.default()
        assert isinstance(str(b), str)


class TestBasalforcingsPicoComprehensive:
    """Comprehensive tests for pico basalforcings class."""

    def test_init(self):
        b = basalforcings.pico()
        assert b is not None

    def test_has_num_basins(self):
        b = basalforcings.pico()
        assert hasattr(b, 'num_basins')
        assert b.num_basins == 0

    def test_has_basin_id(self):
        b = basalforcings.pico()
        assert hasattr(b, 'basin_id')

    def test_has_maxboxcount(self):
        b = basalforcings.pico()
        assert hasattr(b, 'maxboxcount')

    def test_has_overturning_coeff(self):
        b = basalforcings.pico()
        assert hasattr(b, 'overturning_coeff')

    def test_has_gamma_T(self):
        b = basalforcings.pico()
        assert hasattr(b, 'gamma_T')

    def test_has_farocean_temperature(self):
        b = basalforcings.pico()
        assert hasattr(b, 'farocean_temperature')

    def test_has_farocean_salinity(self):
        b = basalforcings.pico()
        assert hasattr(b, 'farocean_salinity')

    def test_has_isplume(self):
        b = basalforcings.pico()
        assert hasattr(b, 'isplume')

    def test_repr(self):
        b = basalforcings.pico()
        r = repr(b)
        assert isinstance(r, str)

    def test_str(self):
        b = basalforcings.pico()
        assert isinstance(str(b), str)


class TestBasalforcingsLinearComprehensive:
    """Comprehensive tests for linear basalforcings class."""

    def test_init(self):
        b = basalforcings.linear()
        assert b is not None

    def test_has_deepwater_melting_rate(self):
        b = basalforcings.linear()
        assert hasattr(b, 'deepwater_melting_rate')

    def test_has_deepwater_elevation(self):
        b = basalforcings.linear()
        assert hasattr(b, 'deepwater_elevation')

    def test_has_upperwater_melting_rate(self):
        b = basalforcings.linear()
        assert hasattr(b, 'upperwater_melting_rate')

    def test_has_upperwater_elevation(self):
        b = basalforcings.linear()
        assert hasattr(b, 'upperwater_elevation')

    def test_has_geothermalflux(self):
        b = basalforcings.linear()
        assert hasattr(b, 'geothermalflux')

    def test_repr(self):
        b = basalforcings.linear()
        r = repr(b)
        assert isinstance(r, str)

    def test_str(self):
        b = basalforcings.linear()
        assert isinstance(str(b), str)


class TestBasalforcingsMismipComprehensive:
    """Comprehensive tests for mismip basalforcings class."""

    def test_init(self):
        b = basalforcings.mismip()
        assert b is not None

    def test_repr(self):
        b = basalforcings.mismip()
        r = repr(b)
        assert isinstance(r, str)

    def test_str(self):
        b = basalforcings.mismip()
        assert isinstance(str(b), str)


class TestBasalforcingsLineararmaComprehensive:
    """Comprehensive tests for lineararma basalforcings class."""

    def test_init(self):
        b = basalforcings.lineararma()
        assert b is not None

    def test_has_num_basins(self):
        b = basalforcings.lineararma()
        assert hasattr(b, 'num_basins')
        assert b.num_basins == 0

    def test_has_basin_id(self):
        b = basalforcings.lineararma()
        assert hasattr(b, 'basin_id')

    def test_has_ar_order(self):
        b = basalforcings.lineararma()
        assert hasattr(b, 'ar_order')

    def test_has_ma_order(self):
        b = basalforcings.lineararma()
        assert hasattr(b, 'ma_order')

    def test_repr(self):
        b = basalforcings.lineararma()
        r = repr(b)
        assert isinstance(r, str)


# plumer class does not exist in current implementation


# sicop class does not exist in current implementation
