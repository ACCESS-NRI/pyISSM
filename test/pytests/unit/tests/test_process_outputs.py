"""
Unit tests for _process_outputs() methods across SMB, hydrology, and other model classes.
These methods are pure string-processing (no ISSM dependency) and can be tested directly.
"""

import pytest
from types import SimpleNamespace

try:
    import pyissm.model.classes.smb as smb
    import pyissm.model.classes.hydrology as hydrology
    # __init__.py shadows module names with class names, so import classes directly
    from pyissm.model.classes.masstransport import masstransport
    from pyissm.model.classes.thermal import thermal
    from pyissm.model.classes.damage import damage
    from pyissm.model.classes.transient import transient
    from pyissm.model.classes.stressbalance import stressbalance
    CLASSES_AVAILABLE = True
except ImportError:
    CLASSES_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not CLASSES_AVAILABLE,
    reason="pyissm model classes not available"
)


def _make_gemb_mesh(ne=10):
    """Minimal mock mesh for smb.gemb() constructor."""
    m = SimpleNamespace()
    m.numberofelements = ne
    return m


# ============== HELPER ==============

def _test_class_process_outputs(cls, expected_default_output):
    """
    Reusable test helper that validates _process_outputs() for a model class.
    Checks:
    - ['default'] expands to the class's default outputs
    - explicit outputs are returned unchanged
    - empty outputs return []
    - return_default_outputs=True returns a tuple
    """
    obj = cls()

    # Default outputs expansion
    obj.requested_outputs = ['default']
    result = obj._process_outputs()
    assert expected_default_output in result, (
        f"{cls.__name__}: expected '{expected_default_output}' in {result}"
    )

    # Custom output passthrough
    obj.requested_outputs = ['CustomOutput']
    result = obj._process_outputs()
    assert result == ['CustomOutput']

    # Empty outputs return empty list
    obj.requested_outputs = []
    result = obj._process_outputs()
    assert result == []

    # Mixed: default + custom
    obj.requested_outputs = ['default', 'CustomOutput']
    result = obj._process_outputs()
    assert expected_default_output in result
    assert 'CustomOutput' in result

    # return_default_outputs=True returns tuple
    obj.requested_outputs = ['default']
    result = obj._process_outputs(return_default_outputs=True)
    assert isinstance(result, tuple)
    assert len(result) == 2
    outputs, defaults = result
    assert expected_default_output in outputs
    assert expected_default_output in defaults


def _make_md_2d():
    """Mock md with a 2D mesh for _process_outputs tests."""
    md = SimpleNamespace()
    md.mesh = SimpleNamespace()
    md.mesh.domain_type = lambda: '2Dhorizontal'
    md.mesh.dimension = lambda: 2
    return md


# ============== SMB CLASSES ==============

class TestSmbDefaultProcessOutputs:
    def test_default_expansion(self):
        _test_class_process_outputs(smb.default, 'SmbMassBalance')

    def test_multiple_outputs(self):
        obj = smb.default()
        obj.requested_outputs = ['default', 'SmbAccumulatedMassBalance']
        result = obj._process_outputs()
        assert 'SmbMassBalance' in result
        assert 'SmbAccumulatedMassBalance' in result

    def test_no_duplicates_for_single_default(self):
        obj = smb.default()
        obj.requested_outputs = ['default']
        result = obj._process_outputs()
        assert result.count('SmbMassBalance') == 1


class TestSmbArmaProcessOutputs:
    def test_default_expansion(self):
        _test_class_process_outputs(smb.arma, 'SmbMassBalance')


class TestSmbComponentsProcessOutputs:
    def test_default_expansion(self):
        _test_class_process_outputs(smb.components, 'SmbMassBalance')


class TestSmbD18opddProcessOutputs:
    def test_default_expansion(self):
        _test_class_process_outputs(smb.d18opdd, 'SmbMassBalance')


class TestSmbGembProcessOutputs:
    def test_default_expansion(self):
        obj = smb.gemb(_make_gemb_mesh())
        obj.requested_outputs = ['default']
        result = obj._process_outputs()
        assert 'SmbMassBalance' in result

    def test_accumulated_in_defaults(self):
        obj = smb.gemb(_make_gemb_mesh())
        obj.requested_outputs = ['default']
        result = obj._process_outputs()
        assert 'SmbAccumulatedMassBalance' in result

    def test_return_default_outputs_has_both(self):
        obj = smb.gemb(_make_gemb_mesh())
        obj.requested_outputs = ['default']
        _, defaults = obj._process_outputs(return_default_outputs=True)
        assert 'SmbMassBalance' in defaults
        assert 'SmbAccumulatedMassBalance' in defaults


class TestSmbGradientsProcessOutputs:
    def test_default_expansion(self):
        _test_class_process_outputs(smb.gradients, 'SmbMassBalance')


class TestSmbGradientscomponentsProcessOutputs:
    def test_default_expansion(self):
        _test_class_process_outputs(smb.gradientscomponents, 'SmbMassBalance')


class TestSmbGradientsElaProcessOutputs:
    def test_default_expansion(self):
        _test_class_process_outputs(smb.gradientsela, 'SmbMassBalance')


class TestSmbHenningProcessOutputs:
    def test_default_expansion(self):
        _test_class_process_outputs(smb.henning, 'SmbMassBalance')


class TestSmbMeltcomponentsProcessOutputs:
    def test_default_expansion(self):
        _test_class_process_outputs(smb.meltcomponents, 'SmbMassBalance')


class TestSmbPddProcessOutputs:
    def test_default_expansion(self):
        _test_class_process_outputs(smb.pdd, 'SmbMassBalance')


class TestSmbPddSicopolisProcessOutputs:
    def test_default_expansion(self):
        _test_class_process_outputs(smb.pddSicopolis, 'SmbMassBalance')


class TestSmbSemicProcessOutputs:
    def test_default_expansion(self):
        _test_class_process_outputs(smb.semic, 'SmbMassBalance')


# ============== HYDROLOGY CLASSES ==============

class TestHydrologyArmapwProcessOutputs:
    def test_explicit_output_passthrough(self):
        obj = hydrology.armapw()
        obj.requested_outputs = ['HydrologyHead']
        result = obj._process_outputs()
        assert result == ['HydrologyHead']

    def test_empty_outputs(self):
        obj = hydrology.armapw()
        obj.requested_outputs = []
        result = obj._process_outputs()
        assert result == []

    def test_default_expands(self):
        obj = hydrology.armapw()
        obj.requested_outputs = ['default']
        result = obj._process_outputs()
        assert isinstance(result, list)
        assert len(result) > 0

    def test_return_default_outputs_true(self):
        obj = hydrology.armapw()
        obj.requested_outputs = ['default']
        result = obj._process_outputs(return_default_outputs=True)
        assert isinstance(result, tuple)


class TestHydrologyDcProcessOutputs:
    def test_explicit_passthrough(self):
        obj = hydrology.dc()
        obj.requested_outputs = ['HydrologyHead']
        assert obj._process_outputs() == ['HydrologyHead']

    def test_default_expands(self):
        obj = hydrology.dc()
        obj.requested_outputs = ['default']
        result = obj._process_outputs()
        assert isinstance(result, list)
        assert len(result) > 0


class TestHydrologyGladsProcessOutputs:
    def test_explicit_passthrough(self):
        obj = hydrology.glads()
        obj.requested_outputs = ['WaterHead']
        assert obj._process_outputs() == ['WaterHead']

    def test_default_expands(self):
        obj = hydrology.glads()
        obj.requested_outputs = ['default']
        result = obj._process_outputs()
        assert isinstance(result, list)


class TestHydrologyPismProcessOutputs:
    def test_explicit_passthrough(self):
        obj = hydrology.pism()
        obj.requested_outputs = ['HydrologicalTillWaterLayer']
        assert obj._process_outputs() == ['HydrologicalTillWaterLayer']

    def test_default_expands(self):
        obj = hydrology.pism()
        obj.requested_outputs = ['default']
        result = obj._process_outputs()
        assert isinstance(result, list)


class TestHydrologyShaktiProcessOutputs:
    def test_explicit_passthrough(self):
        obj = hydrology.shakti()
        obj.requested_outputs = ['HydrologyHead']
        assert obj._process_outputs() == ['HydrologyHead']


class TestHydrologyShreveProcessOutputs:
    def test_explicit_passthrough(self):
        obj = hydrology.shreve()
        obj.requested_outputs = ['HydrologySurface']
        assert obj._process_outputs() == ['HydrologySurface']


class TestHydrologyTwsProcessOutputs:
    def test_explicit_passthrough(self):
        obj = hydrology.tws()
        obj.requested_outputs = ['HydrologyWaterStorage']
        assert obj._process_outputs() == ['HydrologyWaterStorage']


# ============== OTHER MODEL CLASSES ==============

class TestMasstransportProcessOutputs:
    def test_default_expands(self):
        obj = masstransport()
        obj.requested_outputs = ['default']
        result = obj._process_outputs()
        assert isinstance(result, list)
        assert len(result) > 0

    def test_explicit_passthrough(self):
        obj = masstransport()
        obj.requested_outputs = ['Thickness']
        assert obj._process_outputs() == ['Thickness']

    def test_empty(self):
        obj = masstransport()
        obj.requested_outputs = []
        assert obj._process_outputs() == []

    def test_return_default_outputs(self):
        obj = masstransport()
        obj.requested_outputs = ['default']
        result = obj._process_outputs(return_default_outputs=True)
        assert isinstance(result, tuple)


class TestThermalProcessOutputs:
    def test_default_expands(self):
        obj = thermal()
        obj.requested_outputs = ['default']
        result = obj._process_outputs()
        assert isinstance(result, list)
        assert len(result) > 0

    def test_explicit_passthrough(self):
        obj = thermal()
        obj.requested_outputs = ['Temperature']
        assert obj._process_outputs() == ['Temperature']

    def test_return_default_outputs(self):
        obj = thermal()
        obj.requested_outputs = ['default']
        result = obj._process_outputs(return_default_outputs=True)
        assert isinstance(result, tuple)


class TestDamageProcessOutputs:
    def test_default_expands(self):
        obj = damage()
        obj.requested_outputs = ['default']
        result = obj._process_outputs(md=_make_md_2d())
        assert isinstance(result, list)
        assert len(result) > 0

    def test_explicit_passthrough(self):
        obj = damage()
        obj.requested_outputs = ['DamageDbar']
        result = obj._process_outputs(md=_make_md_2d())
        assert 'DamageDbar' in result


class TestTransientProcessOutputs:
    def test_default_expands(self):
        # transient has no built-in default outputs, so expanding 'default' gives []
        obj = transient()
        obj.requested_outputs = ['default']
        result = obj._process_outputs()
        assert isinstance(result, list)

    def test_explicit_passthrough(self):
        obj = transient()
        obj.requested_outputs = ['Thickness']
        assert obj._process_outputs() == ['Thickness']

    def test_return_default_outputs(self):
        obj = transient()
        obj.requested_outputs = ['default']
        result = obj._process_outputs(return_default_outputs=True)
        assert isinstance(result, tuple)


class TestStressbalanceProcessOutputs:
    def test_default_expands(self):
        obj = stressbalance()
        obj.requested_outputs = ['default']
        result = obj._process_outputs(md=_make_md_2d())
        assert isinstance(result, list)
        assert len(result) > 0

    def test_explicit_passthrough(self):
        obj = stressbalance()
        obj.requested_outputs = ['Vel']
        assert obj._process_outputs(md=_make_md_2d()) == ['Vel']

    def test_empty(self):
        obj = stressbalance()
        obj.requested_outputs = []
        assert obj._process_outputs(md=_make_md_2d()) == []

    def test_return_default_outputs(self):
        obj = stressbalance()
        obj.requested_outputs = ['default']
        result = obj._process_outputs(md=_make_md_2d(), return_default_outputs=True)
        assert isinstance(result, tuple)
