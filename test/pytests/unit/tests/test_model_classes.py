"""
Unit tests for pyissm.model.classes module.

Tests cover instantiation, repr, str, and basic functionality of model classes.
"""

import numpy as np
import pytest
from types import SimpleNamespace

try:
    from pyissm.model.classes.constants import constants
    from pyissm.model.classes.debug import debug
    from pyissm.model.classes.private import private
    from pyissm.model.classes.miscellaneous import miscellaneous
    from pyissm.model.classes.mask import mask
    from pyissm.model.classes.geometry import geometry
    from pyissm.model.classes.damage import damage
    from pyissm.model.classes.debris import debris
    from pyissm.model.classes.verbose import verbose
    from pyissm.model.classes.radaroverlay import radaroverlay
    from pyissm.model.classes.outputdefinition import outputdefinition
    from pyissm.model.classes.amr import amr
    from pyissm.model.classes.esa import esa
    from pyissm.model.classes.rotational import rotational
    from pyissm.model.classes.levelset import levelset
    from pyissm.model.classes.steadystate import steadystate
    from pyissm.model.classes.groundingline import groundingline
    from pyissm.model.classes.rifts import rifts
    from pyissm.model.classes.balancethickness import balancethickness
    from pyissm.model.classes.masstransport import masstransport
    from pyissm.model.classes.thermal import thermal
    from pyissm.model.classes.flowequation import flowequation
    from pyissm.model.classes.stressbalance import stressbalance
    from pyissm.model.classes.transient import transient
    from pyissm.model.classes.initialization import initialization
    from pyissm.model.classes.sampling import sampling
    from pyissm.model.classes.surfaceload import surfaceload
    from pyissm.model.classes.lovenumbers import lovenumbers
    from pyissm.model.classes.dependent import dependent
    from pyissm.model.classes.independent import independent
    from pyissm.model.classes.autodiff import autodiff
    from pyissm.model.classes.misfit import misfit
    from pyissm.model.classes.massfluxatgate import massfluxatgate
    from pyissm.model.classes.regionaloutput import regionaloutput
    from pyissm.model.classes.toolkits import toolkits
    from pyissm.model.classes.issmsettings import issmsettings
    from pyissm.model.classes.offlinesolidearthsolution import offlinesolidearthsolution
    CLASSES_AVAILABLE = True
except ImportError:
    CLASSES_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not CLASSES_AVAILABLE,
    reason="Model classes not available"
)


class TestConstants:
    """Tests for constants class."""

    def test_init_defaults(self):
        """Test default initialization."""
        c = constants()
        assert c.g == 9.81
        assert c.omega == 7.292e-5
        assert c.yts == 365.0 * 24.0 * 3600.0
        assert c.referencetemperature == 223.15
        assert c.gravitational_constant == 6.67259e-11

    def test_repr(self):
        """Test string representation."""
        c = constants()
        s = repr(c)
        assert 'gravitational acceleration' in s
        assert 'angular velocity' in s
        assert 'seconds in a year' in s

    def test_str(self):
        """Test short string."""
        c = constants()
        assert 'constants' in str(c).lower()

    def test_init_from_other(self):
        """Test initialization from another object."""
        other = SimpleNamespace(g=10.0, omega=8e-5)
        c = constants(other)
        assert c.g == 10.0
        assert c.omega == 8e-5


class TestDebug:
    """Tests for debug class."""

    def test_init_defaults(self):
        """Test default initialization."""
        d = debug()
        assert d.valgrind == 0
        assert d.gprof == 0
        assert d.profiling == 0

    def test_repr(self):
        """Test string representation."""
        d = debug()
        s = repr(d)
        assert 'valgrind' in s or 'debug' in s.lower()

    def test_str(self):
        """Test short string."""
        d = debug()
        assert 'debug' in str(d).lower()


class TestPrivate:
    """Tests for private class."""

    def test_init_defaults(self):
        """Test default initialization."""
        p = private()
        assert p.isconsistent == 1
        assert p.runtimename == ''
        assert p.solution == ''

    def test_repr(self):
        """Test string representation."""
        p = private()
        s = repr(p)
        assert 'private' in s.lower() or 'isbamg' in s

    def test_str(self):
        """Test short string."""
        p = private()
        assert 'private' in str(p).lower()


class TestMiscellaneous:
    """Tests for miscellaneous class."""

    def test_init_defaults(self):
        """Test default initialization."""
        m = miscellaneous()
        assert m.name == ''
        assert m.notes == ''

    def test_repr(self):
        """Test string representation."""
        m = miscellaneous()
        s = repr(m)
        assert 'name' in s or 'miscellaneous' in s.lower()

    def test_str(self):
        """Test short string."""
        m = miscellaneous()
        assert 'miscellaneous' in str(m).lower()

    def test_set_name(self):
        """Test setting name."""
        m = miscellaneous()
        m.name = 'test_model'
        assert m.name == 'test_model'


class TestMask:
    """Tests for mask class."""

    def test_init_defaults(self):
        """Test default initialization."""
        m = mask()
        assert hasattr(m, 'ocean_levelset')
        assert hasattr(m, 'ice_levelset')

    def test_repr(self):
        """Test string representation."""
        m = mask()
        s = repr(m)
        assert 'levelset' in s or 'mask' in s.lower()

    def test_str(self):
        """Test short string."""
        m = mask()
        assert 'mask' in str(m).lower()


class TestGeometry:
    """Tests for geometry class."""

    def test_init_defaults(self):
        """Test default initialization."""
        g = geometry()
        assert hasattr(g, 'surface')
        assert hasattr(g, 'base')
        assert hasattr(g, 'thickness')
        assert hasattr(g, 'bed')

    def test_repr(self):
        """Test string representation."""
        g = geometry()
        s = repr(g)
        assert 'surface' in s or 'geometry' in s.lower()

    def test_str(self):
        """Test short string."""
        g = geometry()
        assert 'geometry' in str(g).lower()


class TestDamage:
    """Tests for damage class."""

    def test_init_defaults(self):
        """Test default initialization."""
        d = damage()
        assert hasattr(d, 'isdamage')
        assert d.isdamage == 0

    def test_repr(self):
        """Test string representation."""
        d = damage()
        s = repr(d)
        assert 'damage' in s.lower()

    def test_str(self):
        """Test short string."""
        d = damage()
        assert 'damage' in str(d).lower()


class TestDebris:
    """Tests for debris class."""

    def test_init_defaults(self):
        """Test default initialization."""
        d = debris()
        assert hasattr(d, 'spcthickness')
        assert d.stabilization == 2
        assert d.min_thickness == 0.0

    def test_repr(self):
        """Test string representation."""
        d = debris()
        s = repr(d)
        assert 'debris' in s.lower()

    def test_str(self):
        """Test short string."""
        d = debris()
        assert 'debris' in str(d).lower()


class TestVerbose:
    """Tests for verbose class."""

    def test_init_defaults(self):
        """Test default initialization."""
        v = verbose()
        assert hasattr(v, 'solution')
        assert hasattr(v, 'convergence')

    def test_repr(self):
        """Test string representation."""
        v = verbose()
        s = repr(v)
        assert 'verbose' in s.lower() or 'solution' in s

    def test_str(self):
        """Test short string."""
        v = verbose()
        assert 'verbose' in str(v).lower()


class TestRadaroverlay:
    """Tests for radaroverlay class."""

    def test_init_defaults(self):
        """Test default initialization."""
        r = radaroverlay()
        assert hasattr(r, 'pwr')

    def test_repr(self):
        """Test string representation."""
        r = radaroverlay()
        s = repr(r)
        assert 'radar' in s.lower() or 'pwr' in s

    def test_str(self):
        """Test short string."""
        r = radaroverlay()
        assert 'radar' in str(r).lower()


class TestOutputdefinition:
    """Tests for outputdefinition class."""

    def test_init_defaults(self):
        """Test default initialization."""
        o = outputdefinition()
        assert hasattr(o, 'definitions')

    def test_repr(self):
        """Test string representation."""
        o = outputdefinition()
        s = repr(o)
        assert 'output' in s.lower() or 'definition' in s.lower()

    def test_str(self):
        """Test short string."""
        o = outputdefinition()
        assert 'output' in str(o).lower()


class TestAmr:
    """Tests for amr (adaptive mesh refinement) class."""

    def test_init_defaults(self):
        """Test default initialization."""
        a = amr()
        assert hasattr(a, 'hmin')
        assert hasattr(a, 'hmax')

    def test_repr(self):
        """Test string representation."""
        a = amr()
        s = repr(a)
        assert 'amr' in s.lower() or 'mesh' in s.lower() or 'hmin' in s

    def test_str(self):
        """Test short string."""
        a = amr()
        assert 'amr' in str(a).lower()


class TestEsa:
    """Tests for esa (elastic adjustment) class."""

    def test_init_defaults(self):
        """Test default initialization."""
        e = esa()
        assert hasattr(e, 'deltathickness')

    def test_repr(self):
        """Test string representation."""
        e = esa()
        s = repr(e)
        assert 'esa' in s.lower() or 'elastic' in s.lower()

    def test_str(self):
        """Test short string."""
        e = esa()
        assert 'esa' in str(e).lower()


class TestRotational:
    """Tests for rotational class."""

    def test_init_defaults(self):
        """Test default initialization."""
        r = rotational()
        assert hasattr(r, 'equatorialmoi')
        assert hasattr(r, 'polarmoi')

    def test_repr(self):
        """Test string representation."""
        r = rotational()
        s = repr(r)
        assert 'rotational' in s.lower() or 'moi' in s.lower()

    def test_str(self):
        """Test short string."""
        r = rotational()
        assert 'rotational' in str(r).lower()


class TestLevelset:
    """Tests for levelset class."""

    def test_init_defaults(self):
        """Test default initialization."""
        l = levelset()
        assert hasattr(l, 'spclevelset')
        assert hasattr(l, 'reinit_frequency')

    def test_repr(self):
        """Test string representation."""
        l = levelset()
        s = repr(l)
        assert 'levelset' in s.lower()

    def test_str(self):
        """Test short string."""
        l = levelset()
        assert 'levelset' in str(l).lower()


class TestSteadystate:
    """Tests for steadystate class."""

    def test_init_defaults(self):
        """Test default initialization."""
        s = steadystate()
        assert hasattr(s, 'reltol')
        assert hasattr(s, 'maxiter')

    def test_repr(self):
        """Test string representation."""
        ss = steadystate()
        s = repr(ss)
        assert 'steady' in s.lower() or 'reltol' in s

    def test_str(self):
        """Test short string."""
        ss = steadystate()
        assert 'steady' in str(ss).lower()


class TestGroundingline:
    """Tests for groundingline class."""

    def test_init_defaults(self):
        """Test default initialization."""
        g = groundingline()
        assert hasattr(g, 'migration')

    def test_repr(self):
        """Test string representation."""
        g = groundingline()
        s = repr(g)
        assert 'grounding' in s.lower() or 'migration' in s

    def test_str(self):
        """Test short string."""
        g = groundingline()
        assert 'grounding' in str(g).lower()


class TestRifts:
    """Tests for rifts class."""

    def test_init_defaults(self):
        """Test default initialization."""
        r = rifts()
        assert hasattr(r, 'riftstruct')

    def test_repr(self):
        """Test string representation."""
        r = rifts()
        s = repr(r)
        assert 'rift' in s.lower()

    def test_str(self):
        """Test short string."""
        r = rifts()
        assert 'rift' in str(r).lower()


class TestBalancethickness:
    """Tests for balancethickness class."""

    def test_init_defaults(self):
        """Test default initialization."""
        b = balancethickness()
        assert hasattr(b, 'thickening_rate')
        assert hasattr(b, 'stabilization')

    def test_repr(self):
        """Test string representation."""
        b = balancethickness()
        s = repr(b)
        assert 'balance' in s.lower() or 'thickness' in s.lower()

    def test_str(self):
        """Test short string."""
        b = balancethickness()
        assert 'balance' in str(b).lower()


class TestMasstransport:
    """Tests for masstransport class."""

    def test_init_defaults(self):
        """Test default initialization."""
        m = masstransport()
        assert hasattr(m, 'spcthickness')
        assert hasattr(m, 'stabilization')

    def test_repr(self):
        """Test string representation."""
        m = masstransport()
        s = repr(m)
        assert 'mass' in s.lower() or 'transport' in s.lower()

    def test_str(self):
        """Test short string."""
        m = masstransport()
        assert 'mass' in str(m).lower()


class TestThermal:
    """Tests for thermal class."""

    def test_init_defaults(self):
        """Test default initialization."""
        t = thermal()
        assert hasattr(t, 'spctemperature')
        assert hasattr(t, 'stabilization')

    def test_repr(self):
        """Test string representation."""
        t = thermal()
        s = repr(t)
        assert 'thermal' in s.lower() or 'temperature' in s

    def test_str(self):
        """Test short string."""
        t = thermal()
        assert 'thermal' in str(t).lower()


class TestFlowequation:
    """Tests for flowequation class."""

    def test_init_defaults(self):
        """Test default initialization."""
        f = flowequation()
        assert hasattr(f, 'isSIA')
        assert hasattr(f, 'isSSA')
        assert hasattr(f, 'isHO')
        assert hasattr(f, 'isFS')

    def test_repr(self):
        """Test string representation."""
        f = flowequation()
        s = repr(f)
        assert 'flow' in s.lower() or 'equation' in s.lower()

    def test_str(self):
        """Test short string."""
        f = flowequation()
        assert 'flow' in str(f).lower()


class TestStressbalance:
    """Tests for stressbalance class."""

    def test_init_defaults(self):
        """Test default initialization."""
        s = stressbalance()
        assert hasattr(s, 'spcvx')
        assert hasattr(s, 'spcvy')
        assert hasattr(s, 'maxiter')

    def test_repr(self):
        """Test string representation."""
        sb = stressbalance()
        s = repr(sb)
        assert 'stress' in s.lower() or 'balance' in s.lower()

    def test_str(self):
        """Test short string."""
        sb = stressbalance()
        assert 'stress' in str(sb).lower()


class TestTransient:
    """Tests for transient class."""

    def test_init_defaults(self):
        """Test default initialization."""
        t = transient()
        assert hasattr(t, 'ismasstransport')
        assert hasattr(t, 'isstressbalance')
        assert hasattr(t, 'isthermal')

    def test_repr(self):
        """Test string representation."""
        t = transient()
        s = repr(t)
        assert 'transient' in s.lower()

    def test_str(self):
        """Test short string."""
        t = transient()
        assert 'transient' in str(t).lower()


class TestInitialization:
    """Tests for initialization class."""

    def test_init_defaults(self):
        """Test default initialization."""
        i = initialization()
        assert hasattr(i, 'vx')
        assert hasattr(i, 'vy')
        assert hasattr(i, 'temperature')
        assert hasattr(i, 'pressure')

    def test_repr(self):
        """Test string representation."""
        i = initialization()
        s = repr(i)
        assert 'initial' in s.lower()

    def test_str(self):
        """Test short string."""
        i = initialization()
        assert 'initial' in str(i).lower()


class TestSampling:
    """Tests for sampling class."""

    def test_init_defaults(self):
        """Test default initialization."""
        s = sampling()
        assert hasattr(s, 'seed')

    def test_repr(self):
        """Test string representation."""
        sp = sampling()
        s = repr(sp)
        assert 'sampling' in s.lower() or 'seed' in s

    def test_str(self):
        """Test short string."""
        sp = sampling()
        assert 'sampling' in str(sp).lower()


class TestSurfaceload:
    """Tests for surfaceload class."""

    def test_init_defaults(self):
        """Test default initialization."""
        s = surfaceload()
        assert hasattr(s, 'icethicknesschange')

    def test_repr(self):
        """Test string representation."""
        sl = surfaceload()
        s = repr(sl)
        assert 'surface' in s.lower() or 'load' in s.lower()

    def test_str(self):
        """Test short string."""
        sl = surfaceload()
        assert 'surface' in str(sl).lower() or 'load' in str(sl).lower()


class TestLovenumbers:
    """Tests for lovenumbers class."""

    def test_init_defaults(self):
        """Test default initialization."""
        l = lovenumbers()
        assert hasattr(l, 'h')
        assert hasattr(l, 'k')
        assert hasattr(l, 'l')

    def test_repr(self):
        """Test string representation."""
        ln = lovenumbers()
        s = repr(ln)
        assert 'love' in s.lower()

    def test_str(self):
        """Test short string."""
        ln = lovenumbers()
        assert 'love' in str(ln).lower()


class TestDependent:
    """Tests for dependent class."""

    def test_init_defaults(self):
        """Test default initialization."""
        d = dependent()
        assert hasattr(d, 'fos_reverse_index')
        assert d.name == ''
        assert d.index == -1

    def test_repr(self):
        """Test string representation."""
        d = dependent()
        s = repr(d)
        assert 'dependent' in s.lower()

    def test_str(self):
        """Test short string."""
        d = dependent()
        assert 'dependent' in str(d).lower()


class TestIndependent:
    """Tests for independent class."""

    def test_init_defaults(self):
        """Test default initialization."""
        i = independent()
        assert hasattr(i, 'fos_forward_index')

    def test_repr(self):
        """Test string representation."""
        i = independent()
        s = repr(i)
        assert 'independent' in s.lower()

    def test_str(self):
        """Test short string."""
        i = independent()
        assert 'independent' in str(i).lower()


class TestAutodiff:
    """Tests for autodiff class."""

    def test_init_defaults(self):
        """Test default initialization."""
        a = autodiff()
        assert hasattr(a, 'isautodiff')
        assert a.isautodiff == 0

    def test_repr(self):
        """Test string representation."""
        a = autodiff()
        s = repr(a)
        assert 'autodiff' in s.lower() or 'automatic' in s.lower()

    def test_str(self):
        """Test short string."""
        a = autodiff()
        assert 'autodiff' in str(a).lower()


class TestMisfit:
    """Tests for misfit class."""

    def test_init_defaults(self):
        """Test default initialization."""
        m = misfit()
        assert hasattr(m, 'name')
        assert hasattr(m, 'model_string')

    def test_repr(self):
        """Test string representation."""
        m = misfit()
        s = repr(m)
        assert 'misfit' in s.lower()

    def test_str(self):
        """Test short string."""
        m = misfit()
        assert 'misfit' in str(m).lower()


class TestMassfluxatgate:
    """Tests for massfluxatgate class."""

    def test_init_defaults(self):
        """Test default initialization."""
        m = massfluxatgate()
        assert hasattr(m, 'name')

    def test_repr(self):
        """Test string representation."""
        m = massfluxatgate()
        s = repr(m)
        assert 'massflux' in s.lower() or 'gate' in s.lower()

    def test_str(self):
        """Test short string."""
        m = massfluxatgate()
        assert 'massflux' in str(m).lower() or 'gate' in str(m).lower()


class TestRegionaloutput:
    """Tests for regionaloutput class."""

    def test_init_defaults(self):
        """Test default initialization."""
        r = regionaloutput()
        assert hasattr(r, 'name')

    def test_repr(self):
        """Test string representation."""
        r = regionaloutput()
        s = repr(r)
        assert 'regional' in s.lower() or 'output' in s.lower()

    def test_str(self):
        """Test short string."""
        r = regionaloutput()
        assert 'regional' in str(r).lower()


class TestToolkits:
    """Tests for toolkits class."""

    def test_init_defaults(self):
        """Test default initialization."""
        t = toolkits()
        assert hasattr(t, 'DefaultAnalysis')
        assert hasattr(t, 'RecoveryAnalysis')

    def test_repr(self):
        """Test string representation."""
        t = toolkits()
        s = repr(t)
        assert 'toolkit' in s.lower()

    def test_str(self):
        """Test short string."""
        t = toolkits()
        assert 'toolkit' in str(t).lower()


class TestIssmsettings:
    """Tests for issmsettings class."""

    def test_init_defaults(self):
        """Test default initialization."""
        i = issmsettings()
        assert hasattr(i, 'waitonlock')
        assert hasattr(i, 'output_frequency')
        assert i.lowmem == 0

    def test_repr(self):
        """Test string representation."""
        i = issmsettings()
        s = repr(i)
        assert 'settings' in s.lower()

    def test_str(self):
        """Test short string."""
        i = issmsettings()
        assert 'settings' in str(i).lower()


class TestOfflinesolidearthsolution:
    """Tests for offlinesolidearthsolution class."""

    def test_init_defaults(self):
        """Test default initialization."""
        o = offlinesolidearthsolution()
        assert hasattr(o, 'displacementeast')

    def test_repr(self):
        """Test string representation."""
        o = offlinesolidearthsolution()
        s = repr(o)
        assert 'solid' in s.lower() or 'earth' in s.lower()

    def test_str(self):
        """Test short string."""
        o = offlinesolidearthsolution()
        assert 'solid' in str(o).lower() or 'earth' in str(o).lower()
