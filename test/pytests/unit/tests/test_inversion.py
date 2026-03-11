"""
Unit tests for pyissm.model.classes.inversion module.

Tests cover inversion classes: default and m1qn3.
"""

import pytest

try:
    from pyissm.model.classes import inversion
    INVERSION_AVAILABLE = True
except ImportError:
    INVERSION_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not INVERSION_AVAILABLE,
    reason="inversion classes not available"
)


class TestInversionDefault:
    """Tests for default inversion class."""

    def test_init_defaults(self):
        """Test default initialization."""
        i = inversion.default()
        assert i is not None
        assert hasattr(i, 'iscontrol')
        assert hasattr(i, 'incomplete_adjoint')
        assert hasattr(i, 'control_parameters')
        assert hasattr(i, 'nsteps')

    def test_repr(self):
        """Test string representation."""
        i = inversion.default()
        s = repr(i)
        assert isinstance(s, str)
        assert 'inversion' in s.lower()

    def test_str(self):
        """Test short string."""
        i = inversion.default()
        assert isinstance(str(i), str)

    def test_default_values(self):
        """Test default values are set correctly."""
        i = inversion.default()
        assert i.iscontrol == 0
        # incomplete_adjoint can be 0 or 1, just check it exists
        assert i.incomplete_adjoint in [0, 1]


class TestInversionM1qn3:
    """Tests for m1qn3 inversion class."""

    def test_init_defaults(self):
        """Test default initialization."""
        i = inversion.m1qn3()
        assert i is not None
        assert hasattr(i, 'iscontrol')
        assert hasattr(i, 'control_parameters')

    def test_repr(self):
        """Test string representation."""
        i = inversion.m1qn3()
        s = repr(i)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        i = inversion.m1qn3()
        assert isinstance(str(i), str)


class TestInversionSettings:
    """Tests for inversion parameter settings."""

    def test_nsteps_positive(self):
        """Test nsteps is positive."""
        i = inversion.default()
        assert i.nsteps > 0

    def test_control_parameters_is_list(self):
        """Test control_parameters is a list."""
        i = inversion.default()
        assert isinstance(i.control_parameters, list)


class TestInversionInheritance:
    """Tests for inversion class inheritance behavior."""

    def test_default_from_other(self):
        """Test default initialization from another object."""
        i1 = inversion.default()
        i1.iscontrol = 1
        i2 = inversion.default(i1)
        assert i2.iscontrol == i1.iscontrol

    def test_m1qn3_from_other(self):
        """Test m1qn3 initialization from another object."""
        i1 = inversion.m1qn3()
        i1.iscontrol = 1
        i2 = inversion.m1qn3(i1)
        assert i2.iscontrol == i1.iscontrol
