"""
Unit tests for pyissm.model.classes.hydrology module.

Tests cover hydrology classes: dc, glads, shreve, shakti, pism, armapw, tws.
"""

import pytest

try:
    from pyissm.model.classes import hydrology
    HYDROLOGY_AVAILABLE = True
except ImportError:
    HYDROLOGY_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not HYDROLOGY_AVAILABLE,
    reason="hydrology classes not available"
)


class TestHydrologyDc:
    """Tests for DC hydrology class."""

    def test_init_defaults(self):
        """Test default initialization."""
        h = hydrology.dc()
        assert h is not None
        assert hasattr(h, 'water_compressibility')
        assert hasattr(h, 'isefficientlayer')

    def test_repr(self):
        """Test string representation."""
        h = hydrology.dc()
        s = repr(h)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        h = hydrology.dc()
        assert isinstance(str(h), str)


class TestHydrologyGlads:
    """Tests for GlaDS hydrology class."""

    def test_init_defaults(self):
        """Test default initialization."""
        h = hydrology.glads()
        assert h is not None
        assert hasattr(h, 'sheet_conductivity')
        assert hasattr(h, 'ischannels')

    def test_repr(self):
        """Test string representation."""
        h = hydrology.glads()
        s = repr(h)
        assert isinstance(s, str)
        assert 'glads' in s.lower() or 'hydrology' in s.lower()

    def test_str(self):
        """Test short string."""
        h = hydrology.glads()
        assert isinstance(str(h), str)


class TestHydrologyShreve:
    """Tests for Shreve hydrology class."""

    def test_init_defaults(self):
        """Test default initialization."""
        h = hydrology.shreve()
        assert h is not None

    def test_repr(self):
        """Test string representation."""
        h = hydrology.shreve()
        s = repr(h)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        h = hydrology.shreve()
        assert isinstance(str(h), str)


class TestHydrologyShakti:
    """Tests for SHAKTI hydrology class."""

    def test_init_defaults(self):
        """Test default initialization."""
        h = hydrology.shakti()
        assert h is not None

    def test_repr(self):
        """Test string representation."""
        h = hydrology.shakti()
        s = repr(h)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        h = hydrology.shakti()
        assert isinstance(str(h), str)


class TestHydrologyPism:
    """Tests for PISM hydrology class."""

    def test_init_defaults(self):
        """Test default initialization."""
        h = hydrology.pism()
        assert h is not None
        assert hasattr(h, 'drainage_rate')
        assert hasattr(h, 'watercolumn_max')

    def test_repr(self):
        """Test string representation."""
        h = hydrology.pism()
        s = repr(h)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        h = hydrology.pism()
        assert isinstance(str(h), str)


class TestHydrologyArmapw:
    """Tests for ARMA-PW hydrology class."""

    def test_init_defaults(self):
        """Test default initialization."""
        h = hydrology.armapw()
        assert h is not None

    def test_repr(self):
        """Test string representation."""
        h = hydrology.armapw()
        s = repr(h)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        h = hydrology.armapw()
        assert isinstance(str(h), str)


class TestHydrologyTws:
    """Tests for TWS hydrology class."""

    def test_init_defaults(self):
        """Test default initialization."""
        h = hydrology.tws()
        assert h is not None

    def test_repr(self):
        """Test string representation."""
        h = hydrology.tws()
        s = repr(h)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        h = hydrology.tws()
        assert isinstance(str(h), str)


class TestHydrologyInheritance:
    """Tests for hydrology class inheritance behavior."""

    def test_glads_from_other(self):
        """Test GlaDS initialization from another object."""
        h1 = hydrology.glads()
        h1.ischannels = True
        h2 = hydrology.glads(h1)
        assert h2.ischannels == h1.ischannels
