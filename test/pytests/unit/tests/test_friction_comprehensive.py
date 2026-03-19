"""
Comprehensive unit tests for pyissm.model.classes.friction module.

Tests cover all friction parameterization classes.
"""

import pytest

try:
    from pyissm.model.classes import friction
    FRICTION_AVAILABLE = True
except ImportError:
    FRICTION_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not FRICTION_AVAILABLE,
    reason="Friction classes not available"
)


class TestFrictionDefaultComprehensive:
    """Comprehensive tests for default friction class."""

    def test_init(self):
        f = friction.default()
        assert f is not None

    def test_has_coefficient(self):
        f = friction.default()
        assert hasattr(f, 'coefficient')

    def test_has_p(self):
        f = friction.default()
        assert hasattr(f, 'p')

    def test_has_q(self):
        f = friction.default()
        assert hasattr(f, 'q')

    def test_repr(self):
        f = friction.default()
        r = repr(f)
        assert isinstance(r, str)
        assert len(r) > 0

    def test_str(self):
        f = friction.default()
        assert isinstance(str(f), str)


class TestFrictionWeertmanComprehensive:
    """Comprehensive tests for weertman friction class."""

    def test_init(self):
        f = friction.weertman()
        assert f is not None

    def test_has_C(self):
        f = friction.weertman()
        assert hasattr(f, 'C')

    def test_has_m(self):
        f = friction.weertman()
        assert hasattr(f, 'm')

    def test_repr(self):
        f = friction.weertman()
        r = repr(f)
        assert isinstance(r, str)

    def test_str(self):
        f = friction.weertman()
        assert isinstance(str(f), str)


class TestFrictionCoulombComprehensive:
    """Comprehensive tests for coulomb friction class."""

    def test_init(self):
        f = friction.coulomb()
        assert f is not None

    def test_has_coefficient(self):
        f = friction.coulomb()
        assert hasattr(f, 'coefficient')

    def test_has_coefficientcoulomb(self):
        f = friction.coulomb()
        assert hasattr(f, 'coefficientcoulomb')

    def test_has_p(self):
        f = friction.coulomb()
        assert hasattr(f, 'p')

    def test_has_q(self):
        f = friction.coulomb()
        assert hasattr(f, 'q')

    def test_repr(self):
        f = friction.coulomb()
        r = repr(f)
        assert isinstance(r, str)

    def test_str(self):
        f = friction.coulomb()
        assert isinstance(str(f), str)


class TestFrictionSchoofComprehensive:
    """Comprehensive tests for schoof friction class."""

    def test_init(self):
        f = friction.schoof()
        assert f is not None

    def test_has_C(self):
        f = friction.schoof()
        assert hasattr(f, 'C')

    def test_has_m(self):
        f = friction.schoof()
        assert hasattr(f, 'm')

    def test_has_Cmax(self):
        f = friction.schoof()
        assert hasattr(f, 'Cmax')

    def test_has_coupling(self):
        f = friction.schoof()
        assert hasattr(f, 'coupling')
        assert f.coupling == 0

    def test_repr(self):
        f = friction.schoof()
        r = repr(f)
        assert isinstance(r, str)

    def test_str(self):
        f = friction.schoof()
        assert isinstance(str(f), str)


class TestFrictionRegcoulombComprehensive:
    """Comprehensive tests for regcoulomb friction class."""

    def test_init(self):
        f = friction.regcoulomb()
        assert f is not None

    def test_has_C(self):
        f = friction.regcoulomb()
        assert hasattr(f, 'C')

    def test_has_m(self):
        f = friction.regcoulomb()
        assert hasattr(f, 'm')

    def test_repr(self):
        f = friction.regcoulomb()
        r = repr(f)
        assert isinstance(r, str)

    def test_str(self):
        f = friction.regcoulomb()
        assert isinstance(str(f), str)


# budd class does not exist in current implementation


# weertman_temp class does not exist in current implementation


# DNN class does not exist in current implementation


class TestFrictionHydroComprehensive:
    """Comprehensive tests for hydro friction class."""

    def test_init(self):
        f = friction.hydro()
        assert f is not None

    def test_has_coupling(self):
        f = friction.hydro()
        assert hasattr(f, 'coupling')
        assert f.coupling == 0

    def test_has_q(self):
        f = friction.hydro()
        assert hasattr(f, 'q')

    def test_has_C(self):
        f = friction.hydro()
        assert hasattr(f, 'C')

    def test_has_As(self):
        f = friction.hydro()
        assert hasattr(f, 'As')

    def test_repr(self):
        f = friction.hydro()
        r = repr(f)
        assert isinstance(r, str)

    def test_str(self):
        f = friction.hydro()
        assert isinstance(str(f), str)


class TestFrictionWaterlayerComprehensive:
    """Comprehensive tests for waterlayer friction class."""

    def test_init(self):
        f = friction.waterlayer()
        assert f is not None

    def test_has_coefficient(self):
        f = friction.waterlayer()
        assert hasattr(f, 'coefficient')

    def test_has_water_layer(self):
        f = friction.waterlayer()
        assert hasattr(f, 'water_layer')

    def test_has_f(self):
        f = friction.waterlayer()
        assert hasattr(f, 'f')

    def test_repr(self):
        f = friction.waterlayer()
        r = repr(f)
        assert isinstance(r, str)

    def test_str(self):
        f = friction.waterlayer()
        assert isinstance(str(f), str)


# jougeau class does not exist in current implementation
