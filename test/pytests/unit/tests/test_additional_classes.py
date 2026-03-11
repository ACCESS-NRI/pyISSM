"""
Unit tests for additional pyissm.model.classes modules.

Tests cover frontalforcings, dsl, cf, love, qmu, results, and mesh classes.
"""

import pytest

try:
    from pyissm.model.classes import frontalforcings, dsl, love, qmu, results, mesh
    CLASSES_AVAILABLE = True
except ImportError:
    CLASSES_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not CLASSES_AVAILABLE,
    reason="Classes not available"
)


class TestFrontalforcingsDefault:
    """Tests for default frontalforcings class."""

    def test_init_defaults(self):
        """Test default initialization."""
        f = frontalforcings.default()
        assert f is not None
        assert hasattr(f, 'meltingrate')

    def test_repr(self):
        """Test string representation."""
        f = frontalforcings.default()
        s = repr(f)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        f = frontalforcings.default()
        assert isinstance(str(f), str)


class TestFrontalforcingsRignot:
    """Tests for Rignot frontalforcings class."""

    def test_init_defaults(self):
        """Test default initialization."""
        f = frontalforcings.rignot()
        assert f is not None
        # Check actual attributes
        assert hasattr(f, 'subglacial_discharge')

    def test_repr(self):
        """Test string representation."""
        f = frontalforcings.rignot()
        s = repr(f)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        f = frontalforcings.rignot()
        assert isinstance(str(f), str)


class TestFrontalforcingsRignotarma:
    """Tests for Rignot ARMA frontalforcings class."""

    def test_init_defaults(self):
        """Test default initialization."""
        f = frontalforcings.rignotarma()
        assert f is not None
        assert hasattr(f, 'num_basins')

    def test_repr(self):
        """Test string representation."""
        f = frontalforcings.rignotarma()
        s = repr(f)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        f = frontalforcings.rignotarma()
        assert isinstance(str(f), str)


class TestDslDefault:
    """Tests for default DSL class."""

    def test_init_defaults(self):
        """Test default initialization."""
        d = dsl.default()
        assert d is not None
        assert hasattr(d, 'global_average_thermosteric_sea_level')

    def test_repr(self):
        """Test string representation."""
        d = dsl.default()
        s = repr(d)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        d = dsl.default()
        assert isinstance(str(d), str)


class TestDslMme:
    """Tests for MME DSL class."""

    def test_init_defaults(self):
        """Test default initialization."""
        d = dsl.mme()
        assert d is not None
        assert hasattr(d, 'global_average_thermosteric_sea_level')
        assert hasattr(d, 'modelid')

    def test_repr(self):
        """Test string representation."""
        d = dsl.mme()
        s = repr(d)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        d = dsl.mme()
        assert isinstance(str(d), str)


class TestLoveDefault:
    """Tests for default love class."""

    def test_init_defaults(self):
        """Test default initialization."""
        l = love.default()
        assert l is not None
        assert hasattr(l, 'nfreq')
        assert hasattr(l, 'frequencies')
        assert hasattr(l, 'sh_nmin')
        assert hasattr(l, 'sh_nmax')

    def test_repr(self):
        """Test string representation."""
        l = love.default()
        s = repr(l)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        l = love.default()
        assert isinstance(str(l), str)


class TestLoveFourier:
    """Tests for Fourier love class."""

    def test_init_defaults(self):
        """Test default initialization."""
        l = love.fourier()
        assert l is not None
        assert hasattr(l, 'nfreq')
        assert hasattr(l, 'frequencies')

    def test_repr(self):
        """Test string representation."""
        l = love.fourier()
        s = repr(l)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        l = love.fourier()
        assert isinstance(str(l), str)


class TestQmuDefault:
    """Tests for default QMU class."""

    def test_init_defaults(self):
        """Test default initialization."""
        q = qmu.default()
        assert q is not None
        assert hasattr(q, 'isdakota')

    def test_repr(self):
        """Test string representation."""
        q = qmu.default()
        s = repr(q)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        q = qmu.default()
        assert isinstance(str(q), str)


class TestResultsDefault:
    """Tests for default results class."""

    def test_init_defaults(self):
        """Test default initialization."""
        r = results.default()
        assert r is not None

    def test_repr(self):
        """Test string representation."""
        r = results.default()
        s = repr(r)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        r = results.default()
        assert isinstance(str(r), str)


class TestMesh2d:
    """Tests for mesh2d class."""

    def test_init_defaults(self):
        """Test default initialization."""
        m = mesh.mesh2d()
        assert m is not None
        assert hasattr(m, 'x')
        assert hasattr(m, 'y')
        assert hasattr(m, 'elements')

    def test_repr(self):
        """Test string representation."""
        m = mesh.mesh2d()
        s = repr(m)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        m = mesh.mesh2d()
        assert isinstance(str(m), str)


class TestMesh3dprisms:
    """Tests for mesh3dprisms class."""

    def test_init_defaults(self):
        """Test default initialization."""
        m = mesh.mesh3dprisms()
        assert m is not None
        assert hasattr(m, 'x')
        assert hasattr(m, 'y')
        assert hasattr(m, 'z')
        assert hasattr(m, 'elements')

    def test_repr(self):
        """Test string representation."""
        m = mesh.mesh3dprisms()
        s = repr(m)
        assert isinstance(s, str)

    def test_str(self):
        """Test short string."""
        m = mesh.mesh3dprisms()
        assert isinstance(str(m), str)
