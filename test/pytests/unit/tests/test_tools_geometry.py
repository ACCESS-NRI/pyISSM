"""
Unit tests for pyissm.tools.geometry module.

Tests cover:
- Slope computation
- Nowicki ice profile generation

Note: These tests require the ISSM backend to be available.
"""

import numpy as np
import pytest
from types import SimpleNamespace

try:
    from pyissm.tools.geometry import slope, nowicki_profile
    ISSM_AVAILABLE = True
except ImportError:
    ISSM_AVAILABLE = False
    slope = nowicki_profile = None

pytestmark = pytest.mark.skipif(
    not ISSM_AVAILABLE,
    reason="ISSM backend not available"
)


class TestSlope:
    """Tests for the slope function."""

    def test_slope_flat_surface(self, fake_model_2d):
        """Test slope on a flat surface returns zero."""
        # Flat surface: all nodes have same elevation
        fake_model_2d.geometry.surface = np.array([100.0, 100.0, 100.0, 100.0])
        
        sx, sy, s = slope(fake_model_2d)
        
        np.testing.assert_array_almost_equal(sx, np.zeros(2))
        np.testing.assert_array_almost_equal(sy, np.zeros(2))
        np.testing.assert_array_almost_equal(s, np.zeros(2))

    def test_slope_x_direction(self, fake_model_2d):
        """Test slope in x direction only."""
        # Linear gradient in x: increases from left to right
        fake_model_2d.geometry.surface = np.array([0.0, 1.0, 1.0, 0.0])
        
        sx, sy, s = slope(fake_model_2d)
        
        # Both triangles should have positive x-slope
        assert np.all(sx > 0)
        # y-slope should be zero (or very small)
        np.testing.assert_array_almost_equal(sy, np.zeros(2), decimal=10)

    def test_slope_y_direction(self, fake_model_2d):
        """Test slope in y direction only."""
        # Linear gradient in y: increases from bottom to top
        fake_model_2d.geometry.surface = np.array([0.0, 0.0, 1.0, 1.0])
        
        sx, sy, s = slope(fake_model_2d)
        
        # x-slope should be zero
        np.testing.assert_array_almost_equal(sx, np.zeros(2), decimal=10)
        # Both triangles should have positive y-slope
        assert np.all(sy > 0)

    def test_slope_magnitude(self, fake_model_2d):
        """Test slope magnitude computation."""
        # Known gradient for pythagorean triple (3, 4, 5)
        fake_model_2d.geometry.surface = np.array([0.0, 3.0, 7.0, 4.0])
        
        sx, sy, s = slope(fake_model_2d)
        
        # Magnitude should be sqrt(sx^2 + sy^2)
        expected_magnitude = np.sqrt(sx**2 + sy**2)
        np.testing.assert_array_almost_equal(s, expected_magnitude)

    def test_slope_with_custom_field(self, fake_model_2d):
        """Test slope with a custom field (not surface)."""
        custom_field = np.array([10.0, 20.0, 20.0, 10.0])
        
        sx, sy, s = slope(fake_model_2d, field=custom_field)
        
        # Should have non-zero x-slope
        assert np.all(sx > 0)

    def test_slope_returns_correct_shapes(self, fake_model_2d):
        """Test that slope returns arrays of correct shape."""
        sx, sy, s = slope(fake_model_2d)
        
        # Should return one value per element
        num_elements = fake_model_2d.mesh.numberofelements
        assert sx.shape == (num_elements,)
        assert sy.shape == (num_elements,)
        assert s.shape == (num_elements,)

    def test_slope_3d_raises_not_implemented(self, fake_model_3d):
        """Test that 3D mesh raises NotImplementedError."""
        with pytest.raises(NotImplementedError, match="3D meshes not supported"):
            slope(fake_model_3d)


class TestNowickiProfile:
    """Tests for the nowicki_profile function."""

    def test_returns_three_arrays(self):
        """Test that function returns base, thickness, and sea level."""
        x = np.linspace(-100, 100, 201)
        
        b, h, sea = nowicki_profile(x)
        
        assert isinstance(b, np.ndarray)
        assert isinstance(h, np.ndarray)
        assert isinstance(sea, (int, float))

    def test_output_shapes_match_input(self):
        """Test output arrays have same shape as input."""
        x = np.linspace(-100, 100, 201)
        
        b, h, sea = nowicki_profile(x)
        
        assert b.shape == x.shape
        assert h.shape == x.shape

    def test_thickness_positive(self):
        """Test that ice thickness is always positive."""
        x = np.linspace(-100, 100, 201)
        
        b, h, sea = nowicki_profile(x)
        
        # Thickness should be positive everywhere
        assert np.all(h > 0)

    def test_sea_level_positive(self):
        """Test that sea level is positive."""
        x = np.linspace(-100, 100, 201)
        
        b, h, sea = nowicki_profile(x)
        
        assert sea > 0

    def test_grounded_ice_base_at_zero(self):
        """Test that grounded ice (upstream) has base at zero."""
        x = np.linspace(-100, 100, 201)
        
        b, h, sea = nowicki_profile(x)
        
        # First half is grounded (upstream), base should be 0
        mid = len(x) // 2
        np.testing.assert_array_equal(b[:mid], np.zeros(mid))

    def test_floating_ice_base_below_sea(self):
        """Test that floating ice (downstream) has base below sea level."""
        x = np.linspace(-100, 100, 201)
        
        b, h, sea = nowicki_profile(x)
        
        # Second half is floating (downstream)
        mid = len(x) // 2
        # Base should be below sea level for floating ice
        assert np.all(b[mid:] < sea)

    def test_symmetry_not_expected(self):
        """Test profile is not symmetric (different physics upstream/downstream)."""
        x = np.linspace(-100, 100, 200)  # Use even number for equal split
        
        b, h, sea = nowicki_profile(x)
        
        mid = len(x) // 2
        # Upstream and downstream should behave differently
        # (grounded vs floating)
        assert not np.allclose(h[:mid], h[mid:][::-1])

    def test_consistent_results(self):
        """Test that results are consistent across calls."""
        x = np.linspace(-50, 50, 101)
        
        b1, h1, sea1 = nowicki_profile(x)
        b2, h2, sea2 = nowicki_profile(x)
        
        np.testing.assert_array_equal(b1, b2)
        np.testing.assert_array_equal(h1, h2)
        assert sea1 == sea2

    @pytest.mark.parametrize("n_points", [51, 101, 201, 501])
    def test_various_resolutions(self, n_points):
        """Test profile generation at various resolutions."""
        x = np.linspace(-100, 100, n_points)
        
        b, h, sea = nowicki_profile(x)
        
        assert len(b) == n_points
        assert len(h) == n_points
        assert np.all(h > 0)
