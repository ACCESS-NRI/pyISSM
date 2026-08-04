"""
Unit tests for pyissm.model.classes.basalforcings module.

Tests cover basalforcings classes: default, pico, linear, mismip, lineararma.
"""

from types import SimpleNamespace

import numpy as np
import pytest

try:
    from pyissm.model.classes import basalforcings
    BASALFORCINGS_AVAILABLE = True
except ImportError:
    BASALFORCINGS_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not BASALFORCINGS_AVAILABLE,
    reason="basalforcings classes not available"
)


class TestBasalforcingsDefault:
    """Tests for default basalforcings class."""

    def test_init_defaults(self):
        """Test default initialization."""
        b = basalforcings.default()
        assert b is not None

    def test_repr(self):
        """Test string representation."""
        b = basalforcings.default()
        s = repr(b)
        assert isinstance(s, str)
        assert len(s) > 0

    def test_str(self):
        """Test short string."""
        b = basalforcings.default()
        assert isinstance(str(b), str)


class TestBasalforcingsPico:
    """Tests for PICO basalforcings class."""

    def test_init_defaults(self):
        """Test default initialization."""
        b = basalforcings.pico()
        assert b is not None
        assert hasattr(b, 'num_basins')
        assert hasattr(b, 'maxboxcount')
        assert hasattr(b, 'gamma_T')

    def test_repr(self):
        """Test string representation."""
        b = basalforcings.pico()
        s = repr(b)
        assert isinstance(s, str)
        assert 'pico' in s.lower() or 'basal' in s.lower()

    def test_str(self):
        """Test short string."""
        b = basalforcings.pico()
        assert isinstance(str(b), str)


class TestBasalforcingsLinear:
    """Tests for linear basalforcings class."""

    def test_init_defaults(self):
        """Test default initialization."""
        b = basalforcings.linear()
        assert b is not None

    def test_repr(self):
        """Test string representation."""
        b = basalforcings.linear()
        s = repr(b)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        b = basalforcings.linear()
        assert isinstance(str(b), str)


class TestBasalforcingsMismip:
    """Tests for MISMIP basalforcings class."""

    def test_init_defaults(self):
        """Test default initialization."""
        b = basalforcings.mismip()
        assert b is not None

    def test_repr(self):
        """Test string representation."""
        b = basalforcings.mismip()
        s = repr(b)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        b = basalforcings.mismip()
        assert isinstance(str(b), str)


class TestBasalforcingsLineararma:
    """Tests for linear ARMA basalforcings class."""

    def test_init_defaults(self):
        """Test default initialization."""
        b = basalforcings.lineararma()
        assert b is not None
        assert hasattr(b, 'num_basins')

    def test_repr(self):
        """Test string representation."""
        b = basalforcings.lineararma()
        s = repr(b)
        assert isinstance(s, str)


class TestBasalforcingsIsmip6:
    """Tests for ISMIP6 basal forcing fields."""

    def test_extrude_projects_melt_anomaly_as_nodal_timeseries(self):
        mesh3d = SimpleNamespace(
            numberofvertices2d=3,
            numberofelements2d=1,
            numberoflayers=3,
            numberofvertices=9,
            numberofelements=2,
            domain_type=lambda: '3D',
        )
        md = SimpleNamespace(mesh=mesh3d)
        forcing = basalforcings.ismip6()
        forcing.basin_id = np.array([0])
        forcing.geothermalflux = np.array([1.])
        forcing.groundedice_melting_rate = np.ones(3)
        forcing.melt_anomaly = np.array([
            [1., 2.],
            [3., 4.],
            [5., 6.],
            [0., 1.],
        ])

        forcing._extrude(md)

        expected = np.zeros((mesh3d.numberofvertices + 1, 2))
        expected[:mesh3d.numberofvertices2d, :] = np.array([
            [1., 2.],
            [3., 4.],
            [5., 6.],
        ])
        expected[-1, :] = [0., 1.]
        np.testing.assert_array_equal(forcing.melt_anomaly, expected)


class TestBasalforcingsInheritance:
    """Tests for basalforcings class inheritance behavior."""

    def test_pico_from_other(self):
        """Test PICO initialization from another object."""
        b1 = basalforcings.pico()
        b1.num_basins = 5
        b2 = basalforcings.pico(b1)
        assert b2.num_basins == b1.num_basins
