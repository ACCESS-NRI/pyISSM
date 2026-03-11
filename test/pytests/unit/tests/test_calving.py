"""
Unit tests for pyissm.model.classes.calving module.

Tests cover calving classes: default, vonmises, levermann, minthickness, crevassedepth, parameterization.
"""

import pytest

try:
    from pyissm.model.classes import calving
    CALVING_AVAILABLE = True
except ImportError:
    CALVING_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not CALVING_AVAILABLE,
    reason="calving classes not available"
)


class TestCalvingDefault:
    """Tests for default calving class."""

    def test_init_defaults(self):
        """Test default initialization."""
        c = calving.default()
        assert c is not None

    def test_repr(self):
        """Test string representation."""
        c = calving.default()
        s = repr(c)
        assert isinstance(s, str)
        assert len(s) > 0

    def test_str(self):
        """Test short string."""
        c = calving.default()
        assert isinstance(str(c), str)


class TestCalvingVonmises:
    """Tests for von Mises calving class."""

    def test_init_defaults(self):
        """Test default initialization."""
        c = calving.vonmises()
        assert c is not None
        assert hasattr(c, 'stress_threshold_groundedice')
        assert hasattr(c, 'stress_threshold_floatingice')

    def test_repr(self):
        """Test string representation."""
        c = calving.vonmises()
        s = repr(c)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        c = calving.vonmises()
        assert isinstance(str(c), str)


class TestCalvingLevermann:
    """Tests for Levermann calving class."""

    def test_init_defaults(self):
        """Test default initialization."""
        c = calving.levermann()
        assert c is not None
        assert hasattr(c, 'coeff')

    def test_repr(self):
        """Test string representation."""
        c = calving.levermann()
        s = repr(c)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        c = calving.levermann()
        assert isinstance(str(c), str)


class TestCalvingMinthickness:
    """Tests for minimum thickness calving class."""

    def test_init_defaults(self):
        """Test default initialization."""
        c = calving.minthickness()
        assert c is not None
        assert hasattr(c, 'min_thickness')

    def test_repr(self):
        """Test string representation."""
        c = calving.minthickness()
        s = repr(c)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        c = calving.minthickness()
        assert isinstance(str(c), str)


class TestCalvingCrevassedepth:
    """Tests for crevasse depth calving class."""

    def test_init_defaults(self):
        """Test default initialization."""
        c = calving.crevassedepth()
        assert c is not None
        assert hasattr(c, 'crevasse_opening_stress')

    def test_repr(self):
        """Test string representation."""
        c = calving.crevassedepth()
        s = repr(c)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        c = calving.crevassedepth()
        assert isinstance(str(c), str)


class TestCalvingParameterization:
    """Tests for parameterization calving class."""

    def test_init_defaults(self):
        """Test default initialization."""
        c = calving.parameterization()
        assert c is not None
        # Check actual attributes from repr
        assert hasattr(c, 'min_thickness')
        assert hasattr(c, 'use_param')

    def test_repr(self):
        """Test string representation."""
        c = calving.parameterization()
        s = repr(c)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        c = calving.parameterization()
        assert isinstance(str(c), str)


class TestCalvingInheritance:
    """Tests for calving class inheritance behavior."""

    def test_vonmises_from_other(self):
        """Test von Mises initialization from another object."""
        c1 = calving.vonmises()
        c1.stress_threshold_groundedice = 1e6
        c2 = calving.vonmises(c1)
        assert c2.stress_threshold_groundedice == c1.stress_threshold_groundedice
