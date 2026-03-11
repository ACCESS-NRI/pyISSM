"""
Unit tests for pyissm.model.classes.smb module.

Tests cover SMB classes: default, pdd, components, meltcomponents, gradients, gradientsela, arma.
"""

import pytest

try:
    from pyissm.model.classes import smb
    SMB_AVAILABLE = True
except ImportError:
    SMB_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not SMB_AVAILABLE,
    reason="SMB classes not available"
)


class TestSMBDefault:
    """Tests for default SMB class."""

    def test_init_defaults(self):
        """Test default initialization."""
        s = smb.default()
        assert s is not None

    def test_repr(self):
        """Test string representation."""
        s = smb.default()
        r = repr(s)
        assert isinstance(r, str)
        assert len(r) > 0

    def test_str(self):
        """Test short string."""
        s = smb.default()
        assert isinstance(str(s), str)


class TestSMBPdd:
    """Tests for PDD SMB class."""

    def test_init_defaults(self):
        """Test default initialization."""
        s = smb.pdd()
        assert s is not None
        assert hasattr(s, 'desfac')
        assert hasattr(s, 'rlaps')

    def test_repr(self):
        """Test string representation."""
        s = smb.pdd()
        r = repr(s)
        assert isinstance(r, str)

    def test_str(self):
        """Test short string."""
        s = smb.pdd()
        assert isinstance(str(s), str)


class TestSMBComponents:
    """Tests for components SMB class."""

    def test_init_defaults(self):
        """Test default initialization."""
        s = smb.components()
        assert s is not None

    def test_repr(self):
        """Test string representation."""
        s = smb.components()
        r = repr(s)
        assert isinstance(r, str)

    def test_str(self):
        """Test short string."""
        s = smb.components()
        assert isinstance(str(s), str)


class TestSMBMeltcomponents:
    """Tests for meltcomponents SMB class."""

    def test_init_defaults(self):
        """Test default initialization."""
        s = smb.meltcomponents()
        assert s is not None

    def test_repr(self):
        """Test string representation."""
        s = smb.meltcomponents()
        r = repr(s)
        assert isinstance(r, str)

    def test_str(self):
        """Test short string."""
        s = smb.meltcomponents()
        assert isinstance(str(s), str)


class TestSMBGradients:
    """Tests for gradients SMB class."""

    def test_init_defaults(self):
        """Test default initialization."""
        s = smb.gradients()
        assert s is not None
        assert hasattr(s, 'href')
        assert hasattr(s, 'smbref')

    def test_str(self):
        """Test short string."""
        s = smb.gradients()
        assert isinstance(str(s), str)


class TestSMBGradientsela:
    """Tests for gradientsela SMB class."""

    def test_init_defaults(self):
        """Test default initialization."""
        s = smb.gradientsela()
        assert s is not None
        assert hasattr(s, 'ela')

    def test_repr(self):
        """Test string representation."""
        s = smb.gradientsela()
        r = repr(s)
        assert isinstance(r, str)

    def test_str(self):
        """Test short string."""
        s = smb.gradientsela()
        assert isinstance(str(s), str)


class TestSMBArma:
    """Tests for ARMA SMB class."""

    def test_init_defaults(self):
        """Test default initialization."""
        s = smb.arma()
        assert s is not None
        assert hasattr(s, 'num_basins')
        assert hasattr(s, 'ar_order')
        assert hasattr(s, 'ma_order')

    def test_repr(self):
        """Test string representation."""
        s = smb.arma()
        r = repr(s)
        assert isinstance(r, str)

    def test_str(self):
        """Test short string."""
        s = smb.arma()
        assert isinstance(str(s), str)


class TestSMBInheritance:
    """Tests for SMB class inheritance behavior."""

    def test_pdd_from_other(self):
        """Test PDD initialization from another object."""
        s1 = smb.pdd()
        s1.desfac = 0.5
        s2 = smb.pdd(s1)
        assert s2.desfac == s1.desfac
