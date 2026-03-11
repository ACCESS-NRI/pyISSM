"""
Unit tests for pyissm.model.classes.friction module.

Tests cover friction parameterizations: default, weertman, coulomb, schoof, regcoulomb, hydro, waterlayer.
"""

import pytest

try:
    from pyissm.model.classes import friction
    FRICTION_AVAILABLE = True
except ImportError:
    FRICTION_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not FRICTION_AVAILABLE,
    reason="friction classes not available"
)


class TestFrictionDefault:
    """Tests for default friction class."""

    def test_init_defaults(self):
        """Test default initialization."""
        f = friction.default()
        assert f is not None

    def test_repr(self):
        """Test string representation."""
        f = friction.default()
        s = repr(f)
        assert isinstance(s, str)
        assert len(s) > 0

    def test_str(self):
        """Test short string."""
        f = friction.default()
        s = str(f)
        assert isinstance(s, str)


class TestFrictionWeertman:
    """Tests for Weertman friction class."""

    def test_init_defaults(self):
        """Test default initialization."""
        f = friction.weertman()
        assert f is not None
        assert hasattr(f, 'C')
        assert hasattr(f, 'm')

    def test_repr(self):
        """Test string representation."""
        f = friction.weertman()
        s = repr(f)
        assert 'weertman' in s.lower() or 'friction' in s.lower()

    def test_str(self):
        """Test short string."""
        f = friction.weertman()
        s = str(f)
        assert isinstance(s, str)


class TestFrictionCoulomb:
    """Tests for Coulomb friction class."""

    def test_init_defaults(self):
        """Test default initialization."""
        f = friction.coulomb()
        assert f is not None

    def test_repr(self):
        """Test string representation."""
        f = friction.coulomb()
        s = repr(f)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        f = friction.coulomb()
        s = str(f)
        assert isinstance(s, str)


class TestFrictionSchoof:
    """Tests for Schoof friction class."""

    def test_init_defaults(self):
        """Test default initialization."""
        f = friction.schoof()
        assert f is not None
        assert hasattr(f, 'C')
        assert hasattr(f, 'm')

    def test_repr(self):
        """Test string representation."""
        f = friction.schoof()
        s = repr(f)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        f = friction.schoof()
        s = str(f)
        assert isinstance(s, str)


class TestFrictionRegcoulomb:
    """Tests for regularized Coulomb friction class."""

    def test_init_defaults(self):
        """Test default initialization."""
        f = friction.regcoulomb()
        assert f is not None

    def test_repr(self):
        """Test string representation."""
        f = friction.regcoulomb()
        s = repr(f)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        f = friction.regcoulomb()
        s = str(f)
        assert isinstance(s, str)


class TestFrictionHydro:
    """Tests for hydro friction class."""

    def test_init_defaults(self):
        """Test default initialization."""
        f = friction.hydro()
        assert f is not None
        # Check for actual attributes from the repr output
        assert hasattr(f, 'q')
        assert hasattr(f, 'C')
        assert hasattr(f, 'As')
        assert hasattr(f, 'coupling')

    def test_repr(self):
        """Test string representation."""
        f = friction.hydro()
        s = repr(f)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        f = friction.hydro()
        s = str(f)
        assert isinstance(s, str)


class TestFrictionWaterlayer:
    """Tests for waterlayer friction class."""

    def test_init_defaults(self):
        """Test default initialization."""
        f = friction.waterlayer()
        assert f is not None

    def test_repr(self):
        """Test string representation."""
        f = friction.waterlayer()
        s = repr(f)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        f = friction.waterlayer()
        s = str(f)
        assert isinstance(s, str)


class TestFrictionInheritance:
    """Tests for friction class inheritance behavior."""

    def test_weertman_from_other(self):
        """Test Weertman initialization from another object."""
        f1 = friction.weertman()
        f1.C = 1e6
        f2 = friction.weertman(f1)
        assert f2.C == f1.C

    def test_schoof_from_other(self):
        """Test Schoof initialization from another object."""
        f1 = friction.schoof()
        f1.C = 2e6
        f2 = friction.schoof(f1)
        assert f2.C == f1.C
