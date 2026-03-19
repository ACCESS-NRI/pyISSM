"""
Unit tests for pyissm.model.classes - nodalvalue, stochasticforcing, cf, cluster.

Tests cover additional model classes to increase coverage.
"""

import pytest

try:
    from pyissm.model.classes import nodalvalue, stochasticforcing, cf, cluster
    CLASSES_AVAILABLE = True
except ImportError:
    CLASSES_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not CLASSES_AVAILABLE,
    reason="Classes not available"
)


# ============== NODALVALUE TESTS ==============

class TestNodalvalue:
    """Tests for nodalvalue class."""

    def test_init(self):
        """Test initialization."""
        n = nodalvalue.nodalvalue()
        assert n is not None

    def test_has_name(self):
        """Test has name attribute."""
        n = nodalvalue.nodalvalue()
        assert hasattr(n, 'name')

    def test_has_definitionstring(self):
        """Test has definitionstring attribute."""
        n = nodalvalue.nodalvalue()
        assert hasattr(n, 'definitionstring')

    def test_has_model_string(self):
        """Test has model_string attribute."""
        n = nodalvalue.nodalvalue()
        assert hasattr(n, 'model_string')

    def test_has_node(self):
        """Test has node attribute."""
        n = nodalvalue.nodalvalue()
        assert hasattr(n, 'node')

    def test_repr(self):
        """Test string representation."""
        n = nodalvalue.nodalvalue()
        r = repr(n)
        assert isinstance(r, str)

    def test_str(self):
        """Test short string."""
        n = nodalvalue.nodalvalue()
        assert isinstance(str(n), str)


# ============== STOCHASTICFORCING TESTS ==============

class TestStochasticforcing:
    """Tests for stochasticforcing class."""

    def test_init(self):
        """Test initialization."""
        s = stochasticforcing()
        assert s is not None

    def test_has_isstochasticforcing(self):
        """Test has isstochasticforcing attribute."""
        s = stochasticforcing()
        assert hasattr(s, 'isstochasticforcing')
        assert s.isstochasticforcing == 0

    def test_has_fields(self):
        """Test has fields attribute."""
        s = stochasticforcing()
        assert hasattr(s, 'fields')

    def test_has_defaultdimension(self):
        """Test has defaultdimension attribute."""
        s = stochasticforcing()
        assert hasattr(s, 'defaultdimension')
        assert s.defaultdimension == 0

    def test_has_default_id(self):
        """Test has default_id attribute."""
        s = stochasticforcing()
        assert hasattr(s, 'default_id')

    def test_has_covariance(self):
        """Test has covariance attribute."""
        s = stochasticforcing()
        assert hasattr(s, 'covariance')

    def test_has_randomflag(self):
        """Test has randomflag attribute."""
        s = stochasticforcing()
        assert hasattr(s, 'randomflag')
        assert s.randomflag == 1

    def test_repr(self):
        """Test string representation."""
        s = stochasticforcing()
        r = repr(s)
        assert isinstance(r, str)

    def test_str(self):
        """Test short string."""
        s = stochasticforcing()
        assert isinstance(str(s), str)


# ============== CF (COST FUNCTION) TESTS ==============

class TestLevelsetmisfit:
    """Tests for levelsetmisfit class."""

    def test_init(self):
        """Test initialization."""
        c = cf.levelsetmisfit()
        assert c is not None

    def test_has_name(self):
        """Test has name attribute."""
        c = cf.levelsetmisfit()
        assert hasattr(c, 'name')

    def test_has_definitionstring(self):
        """Test has definitionstring attribute."""
        c = cf.levelsetmisfit()
        assert hasattr(c, 'definitionstring')

    def test_has_model_string(self):
        """Test has model_string attribute."""
        c = cf.levelsetmisfit()
        assert hasattr(c, 'model_string')

    def test_has_observation(self):
        """Test has observation attribute."""
        c = cf.levelsetmisfit()
        assert hasattr(c, 'observation')

    def test_has_weights(self):
        """Test has weights attribute."""
        c = cf.levelsetmisfit()
        assert hasattr(c, 'weights')

    def test_repr(self):
        """Test string representation."""
        c = cf.levelsetmisfit()
        r = repr(c)
        assert isinstance(r, str)

    def test_str(self):
        """Test short string."""
        c = cf.levelsetmisfit()
        assert isinstance(str(c), str)


class TestSurfacesquare:
    """Tests for surfacesquare class."""

    def test_init(self):
        """Test initialization."""
        c = cf.surfacesquare()
        assert c is not None

    def test_has_name(self):
        """Test has name attribute."""
        c = cf.surfacesquare()
        assert hasattr(c, 'name')

    def test_has_definitionstring(self):
        """Test has definitionstring attribute."""
        c = cf.surfacesquare()
        assert hasattr(c, 'definitionstring')

    def test_has_datatime(self):
        """Test has datatime attribute."""
        c = cf.surfacesquare()
        assert hasattr(c, 'datatime')

    def test_repr(self):
        """Test string representation."""
        c = cf.surfacesquare()
        r = repr(c)
        assert isinstance(r, str)

    def test_str(self):
        """Test short string."""
        c = cf.surfacesquare()
        assert isinstance(str(c), str)


class TestSurfacesquaretransient:
    """Tests for surfacesquaretransient class."""

    def test_init(self):
        """Test initialization."""
        c = cf.surfacesquaretransient()
        assert c is not None

    def test_has_name(self):
        """Test has name attribute."""
        c = cf.surfacesquaretransient()
        assert hasattr(c, 'name')

    def test_has_observations(self):
        """Test has observations attribute."""
        c = cf.surfacesquaretransient()
        assert hasattr(c, 'observations')

    def test_has_weights(self):
        """Test has weights attribute."""
        c = cf.surfacesquaretransient()
        assert hasattr(c, 'weights')

    def test_repr(self):
        """Test string representation."""
        c = cf.surfacesquaretransient()
        r = repr(c)
        assert isinstance(r, str)

    def test_str(self):
        """Test short string."""
        c = cf.surfacesquaretransient()
        assert isinstance(str(c), str)


class TestSurfacelogvel:
    """Tests for surfacelogvel class."""

    def test_init(self):
        """Test initialization."""
        c = cf.surfacelogvel()
        assert c is not None

    def test_has_name(self):
        """Test has name attribute."""
        c = cf.surfacelogvel()
        assert hasattr(c, 'name')

    def test_has_definitionstring(self):
        """Test has definitionstring attribute."""
        c = cf.surfacelogvel()
        assert hasattr(c, 'definitionstring')

    def test_has_vxobs(self):
        """Test has vxobs attribute."""
        c = cf.surfacelogvel()
        assert hasattr(c, 'vxobs')

    def test_has_vyobs(self):
        """Test has vyobs attribute."""
        c = cf.surfacelogvel()
        assert hasattr(c, 'vyobs')

    def test_repr(self):
        """Test string representation."""
        c = cf.surfacelogvel()
        r = repr(c)
        assert isinstance(r, str)

    def test_str(self):
        """Test short string."""
        c = cf.surfacelogvel()
        assert isinstance(str(c), str)


# ============== CLUSTER TESTS ==============

class TestClusterGeneric:
    """Tests for generic cluster class."""

    def test_init(self):
        """Test initialization."""
        c = cluster.generic()
        assert c is not None

    def test_has_name(self):
        """Test has name attribute."""
        c = cluster.generic()
        assert hasattr(c, 'name')

    def test_has_np(self):
        """Test has np attribute (processors)."""
        c = cluster.generic()
        assert hasattr(c, 'np')
        assert c.np == 1

    def test_has_login(self):
        """Test has login attribute."""
        c = cluster.generic()
        assert hasattr(c, 'login')

    def test_has_port(self):
        """Test has port attribute."""
        c = cluster.generic()
        assert hasattr(c, 'port')
        assert c.port == 0

    def test_has_verbose(self):
        """Test has verbose attribute."""
        c = cluster.generic()
        assert hasattr(c, 'verbose')
        assert c.verbose == 1

    def test_has_codepath(self):
        """Test has codepath attribute."""
        c = cluster.generic()
        assert hasattr(c, 'codepath')

    def test_has_executionpath(self):
        """Test has executionpath attribute."""
        c = cluster.generic()
        assert hasattr(c, 'executionpath')

    def test_has_interactive(self):
        """Test has interactive attribute."""
        c = cluster.generic()
        assert hasattr(c, 'interactive')
        assert c.interactive == 1

    def test_repr(self):
        """Test string representation."""
        c = cluster.generic()
        r = repr(c)
        assert isinstance(r, str)

    def test_str(self):
        """Test short string."""
        c = cluster.generic()
        assert isinstance(str(c), str)


class TestClusterGadi:
    """Tests for gadi cluster class."""

    def test_init(self):
        """Test initialization."""
        c = cluster.gadi()
        assert c is not None

    def test_has_name(self):
        """Test has name attribute."""
        c = cluster.gadi()
        assert hasattr(c, 'name')

    def test_has_np(self):
        """Test has np attribute."""
        c = cluster.gadi()
        assert hasattr(c, 'np')
        assert c.np == 16

    def test_has_time(self):
        """Test has time attribute."""
        c = cluster.gadi()
        assert hasattr(c, 'time')
        assert c.time == 60

    def test_has_memory(self):
        """Test has memory attribute."""
        c = cluster.gadi()
        assert hasattr(c, 'memory')
        assert c.memory == 40

    def test_has_queue(self):
        """Test has queue attribute."""
        c = cluster.gadi()
        assert hasattr(c, 'queue')
        assert c.queue == 'normal'

    def test_has_project(self):
        """Test has project attribute."""
        c = cluster.gadi()
        assert hasattr(c, 'project')

    def test_has_codepath(self):
        """Test has codepath attribute."""
        c = cluster.gadi()
        assert hasattr(c, 'codepath')

    def test_has_executionpath(self):
        """Test has executionpath attribute."""
        c = cluster.gadi()
        assert hasattr(c, 'executionpath')

    def test_has_storage(self):
        """Test has storage attribute."""
        c = cluster.gadi()
        assert hasattr(c, 'storage')

    def test_repr(self):
        """Test string representation."""
        c = cluster.gadi()
        r = repr(c)
        assert isinstance(r, str)

    def test_str(self):
        """Test short string."""
        c = cluster.gadi()
        assert isinstance(str(c), str)

    def test_np_settable(self):
        """Test np can be set."""
        c = cluster.gadi()
        c.np = 48
        assert c.np == 48
