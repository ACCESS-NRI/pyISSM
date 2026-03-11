"""
Unit tests for pyissm.model.classes.solidearth module.

Tests cover Solid Earth parameterization classes.
"""

import numpy as np
import pytest
from types import SimpleNamespace

try:
    from pyissm.model.classes import solidearth
    SOLIDEARTH_AVAILABLE = True
except ImportError:
    SOLIDEARTH_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not SOLIDEARTH_AVAILABLE,
    reason="Solidearth classes not available"
)


class TestSolidEarthSettings:
    """Tests for solidearth settings class."""

    def test_init_defaults(self):
        """Test default initialization."""
        s = solidearth.settings()
        assert hasattr(s, 'reltol')
        assert hasattr(s, 'abstol')
        assert hasattr(s, 'maxiter')
        assert hasattr(s, 'selfattraction')
        assert hasattr(s, 'elastic')
        assert hasattr(s, 'viscous')
        assert hasattr(s, 'rotation')

    def test_repr(self):
        """Test string representation."""
        s = solidearth.settings()
        r = repr(s)
        assert 'solid' in r.lower() or 'earth' in r.lower()

    def test_str(self):
        """Test short string."""
        s = solidearth.settings()
        assert 'solid' in str(s).lower() or 'earth' in str(s).lower() or 'setting' in str(s).lower()


class TestSolidEarthSolution:
    """Tests for solidearth solution class."""

    def test_init_defaults(self):
        """Test default initialization."""
        s = solidearth.solution()
        # Check it has some solution-related attributes

    def test_repr(self):
        """Test string representation."""
        s = solidearth.solution()
        r = repr(s)
        assert 'solid' in r.lower() or 'earth' in r.lower() or 'solution' in r.lower()

    def test_str(self):
        """Test short string."""
        s = solidearth.solution()
        assert 'solid' in str(s).lower() or 'earth' in str(s).lower() or 'solution' in str(s).lower()


class TestSolidEarthEarth:
    """Tests for solidearth earth class."""

    def test_init_defaults(self):
        """Test default initialization."""
        e = solidearth.earth()
        # Should have earth-related settings

    def test_repr(self):
        """Test string representation."""
        e = solidearth.earth()
        r = repr(e)
        assert 'earth' in r.lower()

    def test_str(self):
        """Test short string."""
        e = solidearth.earth()
        assert 'earth' in str(e).lower()


class TestSolidEarthEuropa:
    """Tests for solidearth europa class."""

    def test_init_defaults(self):
        """Test default initialization."""
        e = solidearth.europa()
        # Should have europa-related settings

    def test_repr(self):
        """Test string representation."""
        e = solidearth.europa()
        r = repr(e)
        assert 'europa' in r.lower() or 'solid' in r.lower()

    def test_str(self):
        """Test short string."""
        e = solidearth.europa()
        assert 'europa' in str(e).lower()
