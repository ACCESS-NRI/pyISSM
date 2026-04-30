"""
Unit tests for initialize() methods in pyissm model classes.
These methods set default field values when unset (NaN/zero),
requiring only md.mesh.numberofvertices (and occasionally md.timestepping).
"""

import pytest
import numpy as np
import warnings
from types import SimpleNamespace

try:
    import pyissm.model.classes.smb as smb
    CLASSES_AVAILABLE = True
except ImportError:
    CLASSES_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not CLASSES_AVAILABLE,
    reason="pyissm model classes not available"
)


def _make_md(nv=10):
    """Minimal mock model with mesh and timestepping."""
    md = SimpleNamespace()
    md.mesh = SimpleNamespace()
    md.mesh.numberofvertices = nv
    md.timestepping = SimpleNamespace()
    md.timestepping.time_step = 1.0
    return md


# ============== SMB.DEFAULT ==============

class TestSmbDefaultInitialize:
    def test_returns_self(self):
        obj = smb.default()
        md = _make_md()
        result = obj.initialize(md)
        assert result is obj

    def test_nan_mass_balance_set_to_zeros(self):
        obj = smb.default()
        md = _make_md(nv=5)
        assert np.all(np.isnan(obj.mass_balance))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            obj.initialize(md)
        np.testing.assert_array_equal(obj.mass_balance, np.zeros(5))

    def test_preset_mass_balance_unchanged(self):
        obj = smb.default()
        obj.mass_balance = np.ones(5)
        md = _make_md(nv=5)
        obj.initialize(md)
        np.testing.assert_array_equal(obj.mass_balance, np.ones(5))

    def test_warns_when_mass_balance_set(self):
        obj = smb.default()
        md = _make_md(nv=5)
        with pytest.warns(UserWarning, match='smb.mass_balance not specified'):
            obj.initialize(md)


# ============== SMB.COMPONENTS ==============

class TestSmbComponentsInitialize:
    def test_returns_self(self):
        obj = smb.components()
        md = _make_md()
        result = obj.initialize(md)
        assert result is obj

    def test_nan_accumulation_set_to_zeros(self):
        obj = smb.components()
        md = _make_md(nv=4)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            obj.initialize(md)
        np.testing.assert_array_equal(obj.accumulation, np.zeros(4))

    def test_nan_evaporation_set_to_zeros(self):
        obj = smb.components()
        md = _make_md(nv=4)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            obj.initialize(md)
        np.testing.assert_array_equal(obj.evaporation, np.zeros(4))

    def test_nan_runoff_set_to_zeros(self):
        obj = smb.components()
        md = _make_md(nv=4)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            obj.initialize(md)
        np.testing.assert_array_equal(obj.runoff, np.zeros(4))

    def test_preset_accumulation_unchanged(self):
        obj = smb.components()
        obj.accumulation = np.ones(4)
        md = _make_md(nv=4)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            obj.initialize(md)
        np.testing.assert_array_equal(obj.accumulation, np.ones(4))


# ============== SMB.D18OPDD ==============

class TestSmbD18opddInitialize:
    def test_returns_self(self):
        obj = smb.d18opdd()
        md = _make_md()
        result = obj.initialize(md)
        assert result is obj

    def test_nan_s0p_set_to_zeros(self):
        obj = smb.d18opdd()
        md = _make_md(nv=6)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            obj.initialize(md)
        np.testing.assert_array_equal(obj.s0p, np.zeros(6))

    def test_nan_s0t_set_to_zeros(self):
        obj = smb.d18opdd()
        md = _make_md(nv=6)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            obj.initialize(md)
        np.testing.assert_array_equal(obj.s0t, np.zeros(6))


# ============== SMB.HENNING ==============

class TestSmbHenningInitialize:
    def test_returns_self(self):
        obj = smb.henning()
        md = _make_md()
        result = obj.initialize(md)
        assert result is obj

    def test_nan_smbref_set_to_zeros(self):
        obj = smb.henning()
        md = _make_md(nv=7)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            obj.initialize(md)
        np.testing.assert_array_equal(obj.smbref, np.zeros(7))

    def test_preset_smbref_unchanged(self):
        obj = smb.henning()
        obj.smbref = np.full(7, 3.0)
        md = _make_md(nv=7)
        obj.initialize(md)
        np.testing.assert_array_equal(obj.smbref, np.full(7, 3.0))


# ============== SMB.GEMB (no auto-initialize, just warns) ==============

class TestSmbGembInitialize:
    def test_returns_self(self):
        obj = smb.gemb()
        md = _make_md()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            result = obj.initialize(md)
        assert result is obj

    def test_emits_warning(self):
        obj = smb.gemb()
        md = _make_md()
        with pytest.warns(UserWarning, match='smb.gemb'):
            obj.initialize(md)


# ============== SMB.GRADIENTS ==============

class TestSmbGradientsInitialize:
    def test_returns_self(self):
        obj = smb.gradients()
        md = _make_md()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            result = obj.initialize(md)
        assert result is obj

    def test_emits_warning(self):
        obj = smb.gradients()
        md = _make_md()
        with pytest.warns(UserWarning, match='smb.gradients'):
            obj.initialize(md)


# ============== SMB.GRADIENTSCOMPONENTS ==============

class TestSmbGradientscomponentsInitialize:
    def test_returns_self(self):
        obj = smb.gradientscomponents()
        md = _make_md()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            result = obj.initialize(md)
        assert result is obj

    def test_emits_warning(self):
        obj = smb.gradientscomponents()
        md = _make_md()
        with pytest.warns(UserWarning, match='smb.gradientscomponents'):
            obj.initialize(md)


# ============== SMB.GRADIENTSELA ==============

class TestSmbGradientsElaInitialize:
    def test_returns_self(self):
        obj = smb.gradientsela()
        md = _make_md()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            result = obj.initialize(md)
        assert result is obj

    def test_emits_warning(self):
        obj = smb.gradientsela()
        md = _make_md()
        with pytest.warns(UserWarning, match='smb.gradientsela'):
            obj.initialize(md)


# ============== SMB.ARMA ==============

class TestSmbArmaInitialize:
    def test_returns_self(self):
        obj = smb.arma()
        md = _make_md()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            result = obj.initialize(md)
        assert result is obj

    def test_ar_order_zero_set_to_one(self):
        obj = smb.arma()
        assert obj.ar_order == 0
        md = _make_md()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            obj.initialize(md)
        assert obj.ar_order == 1

    def test_ma_order_zero_set_to_one(self):
        obj = smb.arma()
        assert obj.ma_order == 0
        md = _make_md()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            obj.initialize(md)
        assert obj.ma_order == 1

    def test_arma_timestep_set_from_md(self):
        obj = smb.arma()
        assert obj.arma_timestep == 0
        md = _make_md()
        md.timestepping.time_step = 5.0
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            obj.initialize(md)
        assert obj.arma_timestep == 5.0
