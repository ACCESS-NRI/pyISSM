"""
Unit tests for pyissm.tools.interp module.

Tests cover:
- averaging: mesh-based smoothing function

Note: These tests require the ISSM backend to be available.
"""

import numpy as np
import pytest
from types import SimpleNamespace
from unittest.mock import patch, MagicMock

try:
    from pyissm.tools.interp import averaging
    from pyissm.model.mesh import get_element_areas_volumes
    ISSM_AVAILABLE = True
except ImportError:
    ISSM_AVAILABLE = False
    averaging = None
    get_element_areas_volumes = None

pytestmark = pytest.mark.skipif(
    not ISSM_AVAILABLE,
    reason="ISSM backend not available"
)


def create_simple_2d_mesh():
    """Create a simple 2D mesh model for testing."""
    # Square domain with 4 nodes and 2 triangles
    # 4---3
    # |\  |
    # | \ |
    # |  \|
    # 1---2
    
    mesh = SimpleNamespace()
    mesh.x = np.array([0.0, 1.0, 1.0, 0.0])
    mesh.y = np.array([0.0, 0.0, 1.0, 1.0])
    mesh.elements = np.array([[1, 2, 3], [1, 3, 4]])  # 1-based indexing
    mesh.numberofvertices = 4
    mesh.numberofelements = 2
    mesh.dimension = lambda: 2
    
    md = SimpleNamespace(mesh=mesh)
    return md


def create_single_triangle_mesh():
    """Create a minimal mesh with one triangle."""
    mesh = SimpleNamespace()
    mesh.x = np.array([0.0, 1.0, 0.5])
    mesh.y = np.array([0.0, 0.0, 1.0])
    mesh.elements = np.array([[1, 2, 3]])  # 1-based
    mesh.numberofvertices = 3
    mesh.numberofelements = 1
    mesh.dimension = lambda: 2
    
    md = SimpleNamespace(mesh=mesh)
    return md


class TestAveraging:
    """Tests for the averaging smoothing function."""

    def test_averaging_element_data_zero_iterations(self):
        """Test averaging with element data and 0 iterations."""
        md = create_simple_2d_mesh()
        
        # Element data (one value per element)
        element_data = np.array([1.0, 2.0])
        
        result = averaging(md, element_data, iterations=0)
        
        # Should return nodal values
        assert result.shape == (4, 1)
        assert not np.any(np.isnan(result))

    def test_averaging_nodal_data_zero_iterations(self):
        """Test averaging with nodal data and 0 iterations."""
        md = create_simple_2d_mesh()
        
        # Nodal data (one value per vertex)
        nodal_data = np.array([1.0, 2.0, 3.0, 4.0])
        
        result = averaging(md, nodal_data, iterations=0)
        
        # Should return the same nodal values
        assert result.shape == (4, 1)
        np.testing.assert_array_almost_equal(result.flatten(), nodal_data)

    def test_averaging_smooths_data(self):
        """Test that averaging smooths discontinuous data."""
        md = create_simple_2d_mesh()
        
        # Sharp discontinuity in nodal data
        nodal_data = np.array([10.0, 0.0, 0.0, 0.0])
        
        result_0 = averaging(md, nodal_data, iterations=0)
        result_5 = averaging(md, nodal_data, iterations=5)
        
        # After smoothing, values should be more uniform
        range_0 = np.max(result_0) - np.min(result_0)
        range_5 = np.max(result_5) - np.min(result_5)
        
        assert range_5 < range_0, "Smoothing should reduce data range"

    def test_averaging_invalid_data_length(self):
        """Test that invalid data length raises exception."""
        md = create_simple_2d_mesh()
        
        # Invalid data length (neither nodes nor elements)
        invalid_data = np.array([1.0, 2.0, 3.0])  # 3 values, but 4 nodes and 2 elements
        
        with pytest.raises(Exception):
            averaging(md, invalid_data, iterations=1)

    def test_averaging_preserves_constant(self):
        """Test that averaging preserves constant fields."""
        md = create_simple_2d_mesh()
        
        # Constant nodal data
        constant_value = 5.0
        nodal_data = np.full(4, constant_value)
        
        result = averaging(md, nodal_data, iterations=10)
        
        # All values should remain constant
        np.testing.assert_array_almost_equal(result.flatten(), nodal_data)

    def test_averaging_single_element(self):
        """Test averaging on single element mesh."""
        md = create_single_triangle_mesh()
        
        nodal_data = np.array([1.0, 2.0, 3.0])
        
        result = averaging(md, nodal_data, iterations=1)
        
        assert result.shape == (3, 1)
        assert not np.any(np.isnan(result))

    def test_averaging_element_data_distribution(self):
        """Test that element data is properly distributed to nodes."""
        md = create_simple_2d_mesh()
        
        # Uniform element data
        element_data = np.array([1.0, 1.0])
        
        result = averaging(md, element_data, iterations=0)
        
        # All nodes should have value 1.0
        np.testing.assert_array_almost_equal(
            result.flatten(), 
            np.ones(4),
            decimal=10
        )

    def test_averaging_multiple_iterations(self):
        """Test that multiple iterations produce valid results."""
        md = create_simple_2d_mesh()
        
        nodal_data = np.array([1.0, 2.0, 3.0, 4.0])
        
        for n_iter in [1, 2, 5, 10]:
            result = averaging(md, nodal_data, iterations=n_iter)
            
            assert result.shape == (4, 1)
            assert not np.any(np.isnan(result))
            assert not np.any(np.isinf(result))


class TestAveragingEdgeCases:
    """Edge case tests for averaging function."""

    def test_averaging_zero_data(self):
        """Test averaging with all zeros."""
        md = create_simple_2d_mesh()
        
        zero_data = np.zeros(4)
        
        result = averaging(md, zero_data, iterations=5)
        
        np.testing.assert_array_almost_equal(result.flatten(), zero_data)

    def test_averaging_negative_data(self):
        """Test averaging with negative values."""
        md = create_simple_2d_mesh()
        
        negative_data = np.array([-1.0, -2.0, -3.0, -4.0])
        
        result = averaging(md, negative_data, iterations=3)
        
        assert result.shape == (4, 1)
        assert not np.any(np.isnan(result))

    def test_averaging_large_values(self):
        """Test averaging with large values."""
        md = create_simple_2d_mesh()
        
        large_data = np.array([1e10, 2e10, 3e10, 4e10])
        
        result = averaging(md, large_data, iterations=2)
        
        assert result.shape == (4, 1)
        assert not np.any(np.isnan(result))
        assert not np.any(np.isinf(result))

    def test_averaging_mixed_sign_data(self):
        """Test averaging with mixed positive and negative values."""
        md = create_simple_2d_mesh()
        
        mixed_data = np.array([-10.0, 10.0, -10.0, 10.0])
        
        result = averaging(md, mixed_data, iterations=5)
        
        # Result should be valid (not NaN or Inf)
        assert result.shape == (4, 1)
        assert not np.any(np.isnan(result))
        assert not np.any(np.isinf(result))


class TestAveragingConservation:
    """Tests for conservation properties of averaging."""

    def test_averaging_produces_bounded_results(self):
        """Test that averaging produces results bounded by input range."""
        md = create_simple_2d_mesh()
        
        nodal_data = np.array([1.0, 2.0, 3.0, 4.0])
        
        result = averaging(md, nodal_data, iterations=10)
        
        # Result values should stay within or near the original data range
        assert np.min(result) >= np.min(nodal_data) - 0.5
        assert np.max(result) <= np.max(nodal_data) + 0.5
