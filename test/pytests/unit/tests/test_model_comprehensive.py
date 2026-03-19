"""
Comprehensive unit tests for pyissm.model.Model class.

Tests cover Model initialization, repr, str, and various methods.
"""

import pytest
import copy

try:
    from pyissm.model import Model
    from pyissm.model.classes import class_registry
    MODEL_AVAILABLE = True
except ImportError:
    MODEL_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not MODEL_AVAILABLE,
    reason="Model not available"
)


class TestModelInitializationComprehensive:
    """Comprehensive tests for Model initialization."""

    def test_model_creates_successfully(self):
        """Test that Model() creates without error."""
        md = Model()
        assert md is not None

    def test_model_has_all_expected_components(self):
        """Test Model has all expected sub-components."""
        md = Model()
        expected = [
            'mesh', 'mask', 'geometry', 'constants', 'smb', 'basalforcings',
            'materials', 'damage', 'friction', 'flowequation', 'timestepping',
            'initialization', 'rifts', 'dsl', 'solidearth', 'debug', 'verbose',
            'settings', 'toolkits', 'cluster', 'balancethickness', 'stressbalance',
            'groundingline', 'hydrology', 'debris', 'masstransport', 'thermal',
            'steadystate', 'transient', 'levelset', 'calving', 'frontalforcings',
            'love', 'esa', 'sampling', 'autodiff', 'inversion', 'qmu', 'amr',
            'results', 'outputdefinition', 'radaroverlay', 'miscellaneous',
            'private', 'stochasticforcing'
        ]
        for attr in expected:
            assert hasattr(md, attr), f"Model missing attribute: {attr}"

    def test_repr(self):
        """Test string repr."""
        md = Model()
        r = repr(md)
        assert 'ISSM Model Class' in r
        assert 'mesh' in r
        assert 'geometry' in r

    def test_str(self):
        """Test short string."""
        md = Model()
        assert str(md) == 'ISSM Model Class'


class TestModelClassNames:
    """Tests for model_class_names method."""

    def test_returns_list(self):
        """Test returns a list."""
        md = Model()
        names = md.model_class_names()
        assert isinstance(names, list)

    def test_returns_sorted_list(self):
        """Test returns sorted list."""
        md = Model()
        names = md.model_class_names()
        assert names == sorted(names)

    def test_contains_expected_names(self):
        """Test contains expected component names."""
        md = Model()
        names = md.model_class_names()
        expected = ['mesh', 'mask', 'geometry', 'constants', 'smb']
        for name in expected:
            assert name in names, f"Expected {name} in model_class_names()"


class TestModelCheckMessage:
    """Tests for check_message method."""

    def test_check_message_sets_inconsistent(self, capsys):
        """Test check_message sets isconsistent to False."""
        md = Model()
        md.private.isconsistent = True
        md.check_message("test error")
        assert md.private.isconsistent == False

    def test_check_message_prints_error(self, capsys):
        """Test check_message prints the error."""
        md = Model()
        md.check_message("test consistency error")
        captured = capsys.readouterr()
        assert "Model consistency error" in captured.out
        assert "test consistency error" in captured.out

    def test_check_message_returns_self(self):
        """Test check_message returns self for chaining."""
        md = Model()
        result = md.check_message("some error")
        assert result is md


class TestModelGetState:
    """Tests for __getstate__ method."""

    def test_getstate_returns_dict(self):
        """Test __getstate__ returns a dict."""
        md = Model()
        state = md.__getstate__()
        assert isinstance(state, dict)

    def test_getstate_contains_components(self):
        """Test __getstate__ contains model components."""
        md = Model()
        state = md.__getstate__()
        assert 'mesh' in state
        assert 'geometry' in state
        assert 'mask' in state


class TestModelCopyBehavior:
    """Tests for model copy behavior."""

    def test_deepcopy_creates_independent_model(self):
        """Test deepcopy creates fully independent model."""
        md1 = Model()
        md1.miscellaneous.name = 'original'
        md2 = copy.deepcopy(md1)
        md2.miscellaneous.name = 'copy'
        assert md1.miscellaneous.name == 'original'
        assert md2.miscellaneous.name == 'copy'

    def test_shallow_copy_shares_mutable_objects(self):
        """Test shallow copy behavior."""
        md1 = Model()
        md2 = copy.copy(md1)
        # They are different objects
        assert md1 is not md2


class TestModelMeshAccess:
    """Tests for accessing mesh component."""

    def test_mesh_has_dimension_method(self):
        """Test mesh has dimension method."""
        md = Model()
        assert hasattr(md.mesh, 'dimension')
        assert md.mesh.dimension() == 2

    def test_mesh_has_numberofvertices(self):
        """Test mesh has numberofvertices."""
        md = Model()
        assert hasattr(md.mesh, 'numberofvertices')
        assert md.mesh.numberofvertices == 0

    def test_mesh_has_numberofelements(self):
        """Test mesh has numberofelements."""
        md = Model()
        assert hasattr(md.mesh, 'numberofelements')
        assert md.mesh.numberofelements == 0


class TestModelConstantsAccess:
    """Tests for accessing constants component."""

    def test_constants_has_g(self):
        """Test constants has gravitational acceleration."""
        md = Model()
        assert hasattr(md.constants, 'g')
        assert md.constants.g == 9.81

    def test_constants_has_yts(self):
        """Test constants has years to seconds."""
        md = Model()
        assert hasattr(md.constants, 'yts')


class TestModelTimesteppingAccess:
    """Tests for accessing timestepping component."""

    def test_timestepping_has_start_time(self):
        """Test timestepping has start_time."""
        md = Model()
        assert hasattr(md.timestepping, 'start_time')

    def test_timestepping_has_final_time(self):
        """Test timestepping has final_time."""
        md = Model()
        assert hasattr(md.timestepping, 'final_time')

    def test_timestepping_has_time_step(self):
        """Test timestepping has time_step."""
        md = Model()
        assert hasattr(md.timestepping, 'time_step')


class TestModelSettingsAccess:
    """Tests for accessing settings component."""

    def test_settings_has_results_on_nodes(self):
        """Test settings has results_on_nodes."""
        md = Model()
        assert hasattr(md.settings, 'results_on_nodes')

    def test_settings_has_io_gather(self):
        """Test settings has io_gather."""
        md = Model()
        assert hasattr(md.settings, 'io_gather')


class TestModelTransientAccess:
    """Tests for accessing transient component."""

    def test_transient_has_isthermal(self):
        """Test transient has isthermal."""
        md = Model()
        assert hasattr(md.transient, 'isthermal')

    def test_transient_has_ismasstransport(self):
        """Test transient has ismasstransport."""
        md = Model()
        assert hasattr(md.transient, 'ismasstransport')

    def test_transient_has_isstressbalance(self):
        """Test transient has isstressbalance."""
        md = Model()
        assert hasattr(md.transient, 'isstressbalance')


class TestModelInversionAccess:
    """Tests for accessing inversion component."""

    def test_inversion_has_iscontrol(self):
        """Test inversion has iscontrol."""
        md = Model()
        assert hasattr(md.inversion, 'iscontrol')
        assert md.inversion.iscontrol == 0

    def test_inversion_has_control_parameters(self):
        """Test inversion has control_parameters."""
        md = Model()
        assert hasattr(md.inversion, 'control_parameters')


class TestModelStressbalanceAccess:
    """Tests for accessing stressbalance component."""

    def test_stressbalance_has_spcvx(self):
        """Test stressbalance has spcvx."""
        md = Model()
        assert hasattr(md.stressbalance, 'spcvx')

    def test_stressbalance_has_spcvy(self):
        """Test stressbalance has spcvy."""
        md = Model()
        assert hasattr(md.stressbalance, 'spcvy')


class TestModelMasstransportAccess:
    """Tests for accessing masstransport component."""

    def test_masstransport_has_spcthickness(self):
        """Test masstransport has spcthickness."""
        md = Model()
        assert hasattr(md.masstransport, 'spcthickness')

    def test_masstransport_has_min_thickness(self):
        """Test masstransport has min_thickness."""
        md = Model()
        assert hasattr(md.masstransport, 'min_thickness')


class TestModelThermalAccess:
    """Tests for accessing thermal component."""

    def test_thermal_has_spctemperature(self):
        """Test thermal has spctemperature."""
        md = Model()
        assert hasattr(md.thermal, 'spctemperature')

    def test_thermal_has_stabilization(self):
        """Test thermal has stabilization."""
        md = Model()
        assert hasattr(md.thermal, 'stabilization')


class TestModelMaterialsAccess:
    """Tests for accessing materials component."""

    def test_materials_has_rho_ice(self):
        """Test materials has rho_ice."""
        md = Model()
        assert hasattr(md.materials, 'rho_ice')

    def test_materials_has_rho_water(self):
        """Test materials has rho_water."""
        md = Model()
        assert hasattr(md.materials, 'rho_water')

    def test_materials_has_rheology_B(self):
        """Test materials has rheology_B."""
        md = Model()
        assert hasattr(md.materials, 'rheology_B')

    def test_materials_has_rheology_n(self):
        """Test materials has rheology_n."""
        md = Model()
        assert hasattr(md.materials, 'rheology_n')


class TestModelFrictionAccess:
    """Tests for accessing friction component."""

    def test_friction_has_coefficient(self):
        """Test friction has coefficient."""
        md = Model()
        assert hasattr(md.friction, 'coefficient')

    def test_friction_has_p(self):
        """Test friction has p."""
        md = Model()
        assert hasattr(md.friction, 'p')

    def test_friction_has_q(self):
        """Test friction has q."""
        md = Model()
        assert hasattr(md.friction, 'q')


class TestModelSMBAccess:
    """Tests for accessing SMB component."""

    def test_smb_has_mass_balance(self):
        """Test SMB has mass_balance."""
        md = Model()
        assert hasattr(md.smb, 'mass_balance')


class TestModelCalvingAccess:
    """Tests for accessing calving component."""

    def test_calving_has_calvingrate(self):
        """Test calving has calvingrate."""
        md = Model()
        assert hasattr(md.calving, 'calvingrate')


class TestModelFlowequationAccess:
    """Tests for accessing flowequation component."""

    def test_flowequation_has_isSIA(self):
        """Test flowequation has isSIA."""
        md = Model()
        assert hasattr(md.flowequation, 'isSIA')

    def test_flowequation_has_isSSA(self):
        """Test flowequation has isSSA."""
        md = Model()
        assert hasattr(md.flowequation, 'isSSA')

    def test_flowequation_has_isHO(self):
        """Test flowequation has isHO."""
        md = Model()
        assert hasattr(md.flowequation, 'isHO')

    def test_flowequation_has_isFS(self):
        """Test flowequation has isFS."""
        md = Model()
        assert hasattr(md.flowequation, 'isFS')
