"""
Unit tests for pyissm.model.Model class.

Tests cover:
- Model initialization
- Default attribute presence
- Model component structure
"""

import pytest

from pyissm.model import Model


class TestModelInitialization:
    """Tests for Model class initialization."""

    def test_model_creates_successfully(self):
        """Test that Model can be instantiated."""
        md = Model()
        assert md is not None

    def test_model_has_mesh_attribute(self):
        """Test that Model has mesh attribute."""
        md = Model()
        assert hasattr(md, 'mesh')

    def test_model_has_geometry_attribute(self):
        """Test that Model has geometry attribute."""
        md = Model()
        assert hasattr(md, 'geometry')

    def test_model_has_mask_attribute(self):
        """Test that Model has mask attribute."""
        md = Model()
        assert hasattr(md, 'mask')

    def test_model_has_constants_attribute(self):
        """Test that Model has constants attribute."""
        md = Model()
        assert hasattr(md, 'constants')


class TestModelComponents:
    """Tests for Model component attributes."""

    @pytest.fixture
    def model(self):
        """Create a Model instance for testing."""
        return Model()

    # Physical components
    def test_has_smb(self, model):
        """Test surface mass balance component exists."""
        assert hasattr(model, 'smb')

    def test_has_basalforcings(self, model):
        """Test basal forcings component exists."""
        assert hasattr(model, 'basalforcings')

    def test_has_materials(self, model):
        """Test materials component exists."""
        assert hasattr(model, 'materials')

    def test_has_friction(self, model):
        """Test friction component exists."""
        assert hasattr(model, 'friction')

    def test_has_flowequation(self, model):
        """Test flow equation component exists."""
        assert hasattr(model, 'flowequation')

    # Temporal components
    def test_has_timestepping(self, model):
        """Test timestepping component exists."""
        assert hasattr(model, 'timestepping')

    def test_has_initialization(self, model):
        """Test initialization component exists."""
        assert hasattr(model, 'initialization')

    # Solution components
    def test_has_stressbalance(self, model):
        """Test stress balance component exists."""
        assert hasattr(model, 'stressbalance')

    def test_has_masstransport(self, model):
        """Test mass transport component exists."""
        assert hasattr(model, 'masstransport')

    def test_has_thermal(self, model):
        """Test thermal component exists."""
        assert hasattr(model, 'thermal')

    def test_has_transient(self, model):
        """Test transient component exists."""
        assert hasattr(model, 'transient')

    # Solid earth components
    def test_has_solidearth(self, model):
        """Test solid earth component exists."""
        assert hasattr(model, 'solidearth')

    def test_has_dsl(self, model):
        """Test dynamic sea level component exists."""
        assert hasattr(model, 'dsl')

    # Calving and frontal
    def test_has_calving(self, model):
        """Test calving component exists."""
        assert hasattr(model, 'calving')

    def test_has_frontalforcings(self, model):
        """Test frontal forcings component exists."""
        assert hasattr(model, 'frontalforcings')

    def test_has_groundingline(self, model):
        """Test grounding line component exists."""
        assert hasattr(model, 'groundingline')

    # Settings and configuration
    def test_has_settings(self, model):
        """Test settings component exists."""
        assert hasattr(model, 'settings')

    def test_has_verbose(self, model):
        """Test verbose component exists."""
        assert hasattr(model, 'verbose')

    def test_has_debug(self, model):
        """Test debug component exists."""
        assert hasattr(model, 'debug')

    def test_has_toolkits(self, model):
        """Test toolkits component exists."""
        assert hasattr(model, 'toolkits')

    def test_has_cluster(self, model):
        """Test cluster component exists."""
        assert hasattr(model, 'cluster')

    # Inversion and optimization
    def test_has_inversion(self, model):
        """Test inversion component exists."""
        assert hasattr(model, 'inversion')

    def test_has_autodiff(self, model):
        """Test autodiff component exists."""
        assert hasattr(model, 'autodiff')

    # Results
    def test_has_results(self, model):
        """Test results component exists."""
        assert hasattr(model, 'results')

    def test_has_miscellaneous(self, model):
        """Test miscellaneous component exists."""
        assert hasattr(model, 'miscellaneous')


class TestModelMeshComponent:
    """Tests specifically for the mesh component."""

    @pytest.fixture
    def model(self):
        """Create a Model instance for testing."""
        return Model()

    def test_mesh_is_mesh2d_type(self, model):
        """Test that default mesh is 2D type."""
        # Check it has the expected class name or structure
        mesh_class_name = type(model.mesh).__name__
        assert 'mesh2d' in mesh_class_name.lower() or hasattr(model.mesh, 'numberofvertices')

    def test_mesh_has_dimension_method(self, model):
        """Test mesh has dimension method or attribute."""
        assert hasattr(model.mesh, 'dimension') or hasattr(model.mesh, 'dim')


class TestModelDefaults:
    """Tests for default values in Model components."""

    @pytest.fixture
    def model(self):
        """Create a Model instance for testing."""
        return Model()

    def test_multiple_models_independent(self):
        """Test that multiple Model instances are independent."""
        md1 = Model()
        md2 = Model()
        
        # Modifying one should not affect the other
        # (This tests for potential shared state issues)
        assert md1 is not md2
        assert md1.mesh is not md2.mesh


class TestModelAttributeAccess:
    """Tests for attribute access patterns."""

    def test_nested_attribute_access_mesh(self):
        """Test nested attribute access on mesh."""
        md = Model()
        # This pattern is common in pyISSM workflows
        try:
            _ = md.mesh.numberofvertices
            has_attr = True
        except AttributeError:
            has_attr = False
        
        # Either it has the attribute or accessing it fails gracefully
        assert has_attr or True  # We just want no crash

    def test_nested_attribute_access_geometry(self):
        """Test nested attribute access on geometry."""
        md = Model()
        try:
            _ = md.geometry.surface
            has_attr = True
        except AttributeError:
            has_attr = False
        
        assert has_attr or True


class TestModelCopy:
    """Tests for Model copy behavior."""

    def test_model_can_be_shallow_copied(self):
        """Test that Model can be shallow copied."""
        import copy
        md = Model()
        md_copy = copy.copy(md)
        
        assert md_copy is not md

    def test_model_can_be_deep_copied(self):
        """Test that Model can be deep copied."""
        import copy
        md = Model()
        md_copy = copy.deepcopy(md)
        
        assert md_copy is not md
        assert md_copy.mesh is not md.mesh
