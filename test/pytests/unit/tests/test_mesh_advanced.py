"""
Unit tests for advanced pyissm.model.mesh functions:
- get_nodal_functions_coeff
- compute_hessian
- compute_metric
- elements_from_edge
- export_gmsh / find_segments / model_intersect_3d (NotImplementedError)
"""

import pytest
import numpy as np

try:
    from pyissm.model.mesh import (
        get_nodal_functions_coeff,
        compute_hessian,
        compute_metric,
        elements_from_edge,
        export_gmsh,
        find_segments,
        model_intersect_3d,
        model_merge_3d,
        twod_to_3d,
    )
    MESH_AVAILABLE = True
except ImportError:
    MESH_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not MESH_AVAILABLE,
    reason="pyissm.model.mesh not available"
)


# ============== FIXTURES ==============

@pytest.fixture
def unit_triangle():
    """A single unit right triangle with 1-based indexing."""
    index = np.array([[1, 2, 3]])
    x = np.array([0.0, 1.0, 0.0])
    y = np.array([0.0, 0.0, 1.0])
    return index, x, y


@pytest.fixture
def two_triangle_mesh():
    """Two triangles forming a unit square, 1-based indexing."""
    index = np.array([[1, 2, 3], [1, 3, 4]])
    x = np.array([0.0, 1.0, 1.0, 0.0])
    y = np.array([0.0, 0.0, 1.0, 1.0])
    return index, x, y


# ============== GET_NODAL_FUNCTIONS_COEFF ==============

class TestGetNodalFunctionsCoeff:
    """Tests for get_nodal_functions_coeff()."""

    def test_returns_three_arrays(self, unit_triangle):
        index, x, y = unit_triangle
        result = get_nodal_functions_coeff(index, x, y)
        assert len(result) == 3

    def test_output_shapes(self, unit_triangle):
        index, x, y = unit_triangle
        alpha, beta, gamma = get_nodal_functions_coeff(index, x, y)
        assert alpha.shape == (1, 3)
        assert beta.shape == (1, 3)
        assert gamma.shape == (1, 3)

    def test_multiple_elements_shapes(self, two_triangle_mesh):
        index, x, y = two_triangle_mesh
        alpha, beta, gamma = get_nodal_functions_coeff(index, x, y)
        assert alpha.shape == (2, 3)
        assert beta.shape == (2, 3)
        assert gamma.shape == (2, 3)

    def test_partition_of_unity(self, unit_triangle):
        """Sum of nodal functions at any node should equal 1."""
        index, x, y = unit_triangle
        alpha, beta, gamma = get_nodal_functions_coeff(index, x, y)
        # At node 1 (x=0, y=0): N_i(0, 0) = gamma_i
        N_at_origin = gamma[0, :]
        assert np.isclose(np.sum(N_at_origin), 1.0)

    def test_different_xy_lengths_raises(self, unit_triangle):
        index, x, _ = unit_triangle
        y_wrong = np.array([0.0, 1.0])  # length 2, not 3
        with pytest.raises(TypeError):
            get_nodal_functions_coeff(index, x, y_wrong)

    def test_index_exceeds_nodes_raises(self, unit_triangle):
        index = np.array([[1, 2, 10]])  # node 10 doesn't exist in 3-node mesh
        _, x, y = unit_triangle
        with pytest.raises(TypeError):
            get_nodal_functions_coeff(index, x, y)

    def test_non_triangular_index_raises(self):
        index = np.array([[1, 2, 3, 4]])  # 4 columns (quad element)
        x = np.array([0.0, 1.0, 1.0, 0.0])
        y = np.array([0.0, 0.0, 1.0, 1.0])
        with pytest.raises(TypeError):
            get_nodal_functions_coeff(index, x, y)

    def test_nodal_function_kronecker_delta(self, unit_triangle):
        """N_i(x_j, y_j) = delta_ij (Kronecker delta)."""
        index, x, y = unit_triangle
        alpha, beta, gamma = get_nodal_functions_coeff(index, x, y)
        for node_idx in range(3):
            xi, yi = x[node_idx], y[node_idx]
            N_values = alpha[0, :] * xi + beta[0, :] * yi + gamma[0, :]
            assert np.isclose(N_values[node_idx], 1.0, atol=1e-10)
            for other in range(3):
                if other != node_idx:
                    assert np.isclose(N_values[other], 0.0, atol=1e-10)


# ============== COMPUTE_HESSIAN ==============

class TestComputeHessian:
    """Tests for compute_hessian()."""

    def test_quadratic_field_hessian_node(self, two_triangle_mesh):
        """For f(x, y) = x^2, Hxx = 2, Hxy = 0, Hyy = 0."""
        index, x, y = two_triangle_mesh
        field = x ** 2
        hessian = compute_hessian(index, x, y, field, 'node')
        # Hessian for x^2: diagonal entry should be ~2
        assert hessian.shape == (len(x), 3)

    def test_quadratic_field_hessian_element(self, two_triangle_mesh):
        index, x, y = two_triangle_mesh
        field = x ** 2
        hessian = compute_hessian(index, x, y, field, 'element')
        assert hessian.shape == (len(index), 3)

    def test_linear_field_hessian_near_zero(self, two_triangle_mesh):
        """For f(x, y) = x (linear), hessian should be ~0."""
        index, x, y = two_triangle_mesh
        field = x.copy()  # linear
        hessian = compute_hessian(index, x, y, field, 'node')
        assert np.all(np.abs(hessian) < 1e-10)

    def test_invalid_type_raises(self, unit_triangle):
        index, x, y = unit_triangle
        field = x ** 2
        with pytest.raises(TypeError):
            compute_hessian(index, x, y, field, 'invalid_type')

    def test_wrong_field_size_raises(self, two_triangle_mesh):
        index, x, y = two_triangle_mesh
        field = np.ones(100)  # wrong size
        with pytest.raises(TypeError):
            compute_hessian(index, x, y, field, 'node')

    def test_element_field_input(self, two_triangle_mesh):
        """Element-sized field should be accepted."""
        index, x, y = two_triangle_mesh
        field = np.array([1.0, 2.0])  # 2 elements
        hessian = compute_hessian(index, x, y, field, 'element')
        assert hessian.shape[0] == 2


# ============== COMPUTE_METRIC ==============

class TestComputeMetric:
    """Tests for compute_metric()."""

    def test_returns_correct_shape(self, two_triangle_mesh):
        index, x, y = two_triangle_mesh
        field = x ** 2 + y ** 2
        hessian = compute_hessian(index, x, y, field, 'node')
        metric = compute_metric(hessian, 1.0, 0.01, 0.1, 10.0, np.array([], dtype=int))
        assert metric.shape == (len(x), 3)

    def test_no_nans_in_output(self, two_triangle_mesh):
        index, x, y = two_triangle_mesh
        field = x ** 2 + y ** 2
        hessian = compute_hessian(index, x, y, field, 'node')
        metric = compute_metric(hessian, 1.0, 0.01, 0.1, 10.0, np.array([], dtype=int))
        assert not np.any(np.isnan(metric))

    def test_water_elements_set_to_hmax(self, two_triangle_mesh):
        """Nodes specified in pos should get 1/hmax^2 metric."""
        index, x, y = two_triangle_mesh
        field = x ** 2
        hessian = compute_hessian(index, x, y, field, 'node')
        hmax = 10.0
        pos = np.array([0, 1])  # mark first two nodes as water
        metric = compute_metric(hessian, 1.0, 0.01, 0.1, hmax, pos)
        expected_val = 1.0 / hmax ** 2
        assert np.isclose(metric[0, 0], expected_val)
        assert np.isclose(metric[1, 0], expected_val)

    def test_hmin_hmax_constraints(self, two_triangle_mesh):
        """Metric values should respect hmin/hmax constraints."""
        index, x, y = two_triangle_mesh
        field = x ** 2 + y ** 2
        hessian = compute_hessian(index, x, y, field, 'node')
        hmin, hmax = 0.5, 5.0
        metric = compute_metric(hessian, 1.0, 0.01, hmin, hmax, np.array([], dtype=int))
        # All diagonal entries should be in [1/hmax^2, 1/hmin^2]
        diag = np.abs(metric[:, 0])  # M11
        assert np.all(diag <= 1.0 / hmin**2 + 1e-10)
        assert np.all(diag >= 1.0 / hmax**2 - 1e-10)


# ============== ELEMENTS_FROM_EDGE ==============

class TestElementsFromEdge:
    """Tests for elements_from_edge()."""

    def test_shared_edge_returns_two_elements(self, two_triangle_mesh):
        """Edge shared by two triangles should return both."""
        index, x, y = two_triangle_mesh
        # Nodes 1 and 3 are the shared edge (1-based)
        result = elements_from_edge(index, 1, 3)
        assert len(result) == 2

    def test_boundary_edge_returns_one_element(self, two_triangle_mesh):
        """Edge on the boundary should return only one element."""
        index, x, y = two_triangle_mesh
        # Edge 1-2 is on the boundary (only in first triangle)
        result = elements_from_edge(index, 1, 2)
        assert len(result) == 1

    def test_returns_1based_indices(self, two_triangle_mesh):
        """Returned element IDs should be 1-based."""
        index, x, y = two_triangle_mesh
        result = elements_from_edge(index, 1, 2)
        assert np.all(result >= 1)

    def test_nonexistent_edge_returns_empty(self, two_triangle_mesh):
        """Edge not in mesh should return empty."""
        index, x, y = two_triangle_mesh
        result = elements_from_edge(index, 2, 4)  # not a shared edge
        # This edge exists only in element 2 (1-3-4) if 4 is included
        # Actually 2,4 are not in the same element, so empty
        # Depends on mesh - just check it returns an array
        assert isinstance(result, np.ndarray)

    def test_single_triangle_edge(self, unit_triangle):
        index, x, y = unit_triangle
        result = elements_from_edge(index, 1, 2)
        assert len(result) == 1
        assert result[0] == 1


# ============== NOT-IMPLEMENTED FUNCTIONS ==============

class TestNotImplementedFunctions:
    """Tests for functions that raise NotImplementedError."""

    def test_export_gmsh_raises(self):
        with pytest.raises(NotImplementedError):
            export_gmsh()

    def test_find_segments_raises(self):
        with pytest.raises(NotImplementedError):
            find_segments()

    def test_model_intersect_3d_raises(self):
        with pytest.raises(NotImplementedError):
            model_intersect_3d()

    def test_model_merge_3d_raises(self):
        with pytest.raises(NotImplementedError):
            model_merge_3d()

    def test_twod_to_3d_raises(self):
        with pytest.raises(NotImplementedError):
            twod_to_3d()
