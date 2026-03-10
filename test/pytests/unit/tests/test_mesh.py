"""
Unit tests for pyissm.model.mesh module.

Tests cover:
- Mesh creation (get_mesh)
- Mesh processing (process_mesh)
- Node type identification (find_node_types)
"""

import numpy as np
import pytest
from types import SimpleNamespace
import matplotlib.tri as tri

from pyissm.model.mesh import get_mesh, process_mesh, find_node_types


class TestGetMesh:
    """Tests for the get_mesh function."""

    def test_creates_triangulation(self):
        """Test that get_mesh creates a matplotlib Triangulation."""
        x = np.array([0.0, 1.0, 0.5])
        y = np.array([0.0, 0.0, 1.0])
        elements = np.array([[0, 1, 2]])
        
        mesh = get_mesh(x, y, elements)
        
        assert isinstance(mesh, tri.Triangulation)

    def test_handles_1based_indexing(self):
        """Test that get_mesh converts 1-based to 0-based indexing."""
        x = np.array([0.0, 1.0, 0.5])
        y = np.array([0.0, 0.0, 1.0])
        elements_1based = np.array([[1, 2, 3]])  # 1-based
        
        mesh = get_mesh(x, y, elements_1based)
        
        # Should not raise and should have valid triangles
        assert mesh.triangles.min() == 0
        assert mesh.triangles.max() == 2

    def test_handles_0based_indexing(self):
        """Test that get_mesh works with 0-based indexing."""
        x = np.array([0.0, 1.0, 0.5])
        y = np.array([0.0, 0.0, 1.0])
        elements_0based = np.array([[0, 1, 2]])
        
        mesh = get_mesh(x, y, elements_0based)
        
        assert mesh.triangles.min() == 0

    def test_multiple_triangles(self):
        """Test mesh with multiple triangles."""
        x = np.array([0.0, 1.0, 1.0, 0.0])
        y = np.array([0.0, 0.0, 1.0, 1.0])
        elements = np.array([[1, 2, 3], [1, 3, 4]])
        
        mesh = get_mesh(x, y, elements)
        
        assert len(mesh.triangles) == 2

    def test_mesh_coordinates_preserved(self):
        """Test that mesh coordinates are preserved."""
        x = np.array([0.0, 1.0, 0.5])
        y = np.array([0.0, 0.0, 1.0])
        elements = np.array([[0, 1, 2]])
        
        mesh = get_mesh(x, y, elements)
        
        np.testing.assert_array_equal(mesh.x, x)
        np.testing.assert_array_equal(mesh.y, y)


class TestProcessMesh:
    """Tests for the process_mesh function."""

    def test_processes_2d_mesh(self, fake_model_2d):
        """Test processing a 2D mesh."""
        mesh, mesh_x, mesh_y, mesh_elements, is3d = process_mesh(fake_model_2d)
        
        assert isinstance(mesh, tri.Triangulation)
        assert is3d is False
        assert len(mesh_x) == fake_model_2d.mesh.numberofvertices

    def test_returns_0based_elements(self, fake_model_2d):
        """Test that returned elements are 0-based."""
        mesh, mesh_x, mesh_y, mesh_elements, is3d = process_mesh(fake_model_2d)
        
        assert mesh_elements.min() == 0

    def test_processes_3d_mesh(self, fake_model_3d):
        """Test processing a 3D mesh extracts 2D layer."""
        import warnings
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            mesh, mesh_x, mesh_y, mesh_elements, is3d = process_mesh(fake_model_3d)
            
            # Should warn about 3D model
            assert any("3D model" in str(warning.message) for warning in w)
        
        assert is3d is True
        # Should use 2D coordinates
        assert len(mesh_x) == fake_model_3d.mesh.numberofvertices2d

    def test_mesh_x_y_match_elements(self, fake_model_2d):
        """Test that coordinates are consistent with elements."""
        mesh, mesh_x, mesh_y, mesh_elements, is3d = process_mesh(fake_model_2d)
        
        # All element indices should be valid
        assert mesh_elements.max() < len(mesh_x)
        assert mesh_elements.max() < len(mesh_y)


class TestFindNodeTypes:
    """Tests for the find_node_types function."""

    def test_returns_dict(self, fake_model_2d):
        """Test that function returns a dictionary."""
        ice_ls = fake_model_2d.mask.ice_levelset
        ocean_ls = fake_model_2d.mask.ocean_levelset
        
        result = find_node_types(fake_model_2d, ice_ls, ocean_ls)
        
        assert isinstance(result, dict)

    def test_returns_all_keys(self, fake_model_2d):
        """Test that all expected keys are present."""
        ice_ls = fake_model_2d.mask.ice_levelset
        ocean_ls = fake_model_2d.mask.ocean_levelset
        
        result = find_node_types(fake_model_2d, ice_ls, ocean_ls)
        
        expected_keys = ['ice_nodes', 'ice_front_nodes', 'ocean_nodes', 
                        'floating_ice_nodes', 'grounded_ice_nodes']
        for key in expected_keys:
            assert key in result

    def test_all_ice_nodes(self, fake_model_2d):
        """Test identification when all nodes have ice."""
        # All ice (negative level set)
        ice_ls = np.array([-1.0, -1.0, -1.0, -1.0])
        ocean_ls = np.array([1.0, 1.0, 1.0, 1.0])  # No ocean
        
        result = find_node_types(fake_model_2d, ice_ls, ocean_ls)
        
        assert np.all(result['ice_nodes'])
        assert not np.any(result['ice_front_nodes'])

    def test_no_ice_nodes(self, fake_model_2d):
        """Test identification when no nodes have ice."""
        # No ice (positive level set)
        ice_ls = np.array([1.0, 1.0, 1.0, 1.0])
        ocean_ls = np.array([-1.0, -1.0, -1.0, -1.0])  # All ocean
        
        result = find_node_types(fake_model_2d, ice_ls, ocean_ls)
        
        assert not np.any(result['ice_nodes'])
        assert np.all(result['ocean_nodes'])

    def test_ice_front_at_zero(self, fake_model_2d):
        """Test that ice front is identified at zero level set."""
        ice_ls = np.array([-1.0, 0.0, -1.0, 0.0])  # Two nodes at front
        ocean_ls = np.array([1.0, 1.0, 1.0, 1.0])
        
        result = find_node_types(fake_model_2d, ice_ls, ocean_ls)
        
        assert result['ice_front_nodes'][1] is True or result['ice_front_nodes'][1] == True
        assert result['ice_front_nodes'][3] is True or result['ice_front_nodes'][3] == True

    def test_floating_ice_detection(self, fake_model_2d):
        """Test detection of floating ice (ice & ocean overlap)."""
        ice_ls = np.array([-1.0, -1.0, -1.0, -1.0])   # All ice
        ocean_ls = np.array([-1.0, -1.0, 1.0, 1.0])   # Half ocean
        
        result = find_node_types(fake_model_2d, ice_ls, ocean_ls)
        
        # Floating where both ice and ocean
        np.testing.assert_array_equal(
            result['floating_ice_nodes'], 
            [True, True, False, False]
        )

    def test_grounded_ice_detection(self, fake_model_2d):
        """Test detection of grounded ice (ice & no ocean)."""
        ice_ls = np.array([-1.0, -1.0, -1.0, -1.0])   # All ice
        ocean_ls = np.array([-1.0, -1.0, 1.0, 1.0])   # Half ocean
        
        result = find_node_types(fake_model_2d, ice_ls, ocean_ls)
        
        # Grounded where ice but no ocean
        np.testing.assert_array_equal(
            result['grounded_ice_nodes'], 
            [False, False, True, True]
        )

    def test_floating_grounded_mutually_exclusive(self, fake_model_2d):
        """Test that floating and grounded are mutually exclusive for ice."""
        ice_ls = np.array([-1.0, -1.0, -1.0, -1.0])
        ocean_ls = np.array([-1.0, 1.0, -1.0, 1.0])
        
        result = find_node_types(fake_model_2d, ice_ls, ocean_ls)
        
        # No node should be both floating and grounded
        overlap = result['floating_ice_nodes'] & result['grounded_ice_nodes']
        assert not np.any(overlap)

    def test_result_arrays_correct_length(self, fake_model_2d):
        """Test that result arrays have correct length."""
        ice_ls = fake_model_2d.mask.ice_levelset
        ocean_ls = fake_model_2d.mask.ocean_levelset
        
        result = find_node_types(fake_model_2d, ice_ls, ocean_ls)
        
        for key in result:
            assert len(result[key]) == fake_model_2d.mesh.numberofvertices


class TestMeshEdgeCases:
    """Edge case tests for mesh functions."""

    def test_single_triangle_mesh(self):
        """Test mesh functions work with a single triangle."""
        mesh = SimpleNamespace()
        mesh.x = np.array([0.0, 1.0, 0.5])
        mesh.y = np.array([0.0, 0.0, 1.0])
        mesh.elements = np.array([[1, 2, 3]])
        mesh.numberofvertices = 3
        mesh.numberofelements = 1
        mesh.dimension = lambda: 2
        
        md = SimpleNamespace(mesh=mesh)
        md.geometry = SimpleNamespace()
        md.mask = SimpleNamespace()
        
        # Should work without error
        tri_mesh, x, y, elements, is3d = process_mesh(md)
        assert len(tri_mesh.triangles) == 1

    def test_large_mesh(self):
        """Test mesh functions work with larger meshes."""
        # Create a grid of points
        n = 10
        x = np.linspace(0, 1, n)
        y = np.linspace(0, 1, n)
        xx, yy = np.meshgrid(x, y)
        points_x = xx.flatten()
        points_y = yy.flatten()
        
        # Create simple triangulation
        tri_obj = tri.Triangulation(points_x, points_y)
        elements = tri_obj.triangles + 1  # Convert to 1-based
        
        mesh = get_mesh(points_x, points_y, elements)
        
        assert len(mesh.x) == n * n
        assert len(mesh.triangles) > 0
