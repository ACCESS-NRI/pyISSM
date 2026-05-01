"""
Unit tests for pyissm.model.io - _collapse_solution_to_step and _expand_step_to_solution.
These are pure Python functions (no ISSM/netCDF4 dependency at runtime).
"""

import pytest
import numpy as np

try:
    from pyissm.model.io import _collapse_solution_to_step, _expand_step_to_solution
    from pyissm.model.classes.results import solution, solutionstep
    IO_AVAILABLE = True
except ImportError:
    IO_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not IO_AVAILABLE,
    reason="pyissm.model.io not available (likely missing netCDF4)"
)


def _make_step(time, step_num, **kwargs):
    """Helper: create a solutionstep with given attributes."""
    s = solutionstep()
    s.time = float(time)
    s.step = int(step_num)
    for k, v in kwargs.items():
        setattr(s, k, v)
    return s


def _make_solution(steps):
    """Helper: wrap a list of solutionstep in a solution object."""
    sol = solution()
    sol.steps = steps
    return sol


# ============== _COLLAPSE_SOLUTION_TO_STEP ==============

class TestCollapseSolutionToStep:
    """Tests for _collapse_solution_to_step()."""

    def test_empty_solution_returns_empty_step(self):
        sol = solution()
        sol.steps = []
        result = _collapse_solution_to_step(sol)
        assert isinstance(result, solutionstep)

    def test_single_step_time_as_array(self):
        steps = [_make_step(1.0, 1, Thickness=np.array([1.0, 2.0, 3.0]))]
        sol = _make_solution(steps)
        result = _collapse_solution_to_step(sol)
        assert hasattr(result, 'time')
        assert isinstance(result.time, np.ndarray)
        assert result.time[0] == pytest.approx(1.0)

    def test_multiple_steps_time_array(self):
        steps = [
            _make_step(1.0, 1, Vel=np.array([1.0, 2.0])),
            _make_step(2.0, 2, Vel=np.array([3.0, 4.0])),
            _make_step(3.0, 3, Vel=np.array([5.0, 6.0])),
        ]
        sol = _make_solution(steps)
        result = _collapse_solution_to_step(sol)
        np.testing.assert_array_equal(result.time, [1.0, 2.0, 3.0])
        np.testing.assert_array_equal(result.step, [1, 2, 3])

    def test_arrays_stacked_along_time_axis(self):
        steps = [
            _make_step(1.0, 1, Thickness=np.array([1.0, 2.0])),
            _make_step(2.0, 2, Thickness=np.array([3.0, 4.0])),
        ]
        sol = _make_solution(steps)
        result = _collapse_solution_to_step(sol)
        assert hasattr(result, 'Thickness')
        # stacked shape: (2, 2) -> squeeze -> (2, 2) or (2,)
        assert result.Thickness.ndim >= 1

    def test_step_values_collected(self):
        steps = [_make_step(1.0, 1), _make_step(2.0, 2)]
        sol = _make_solution(steps)
        result = _collapse_solution_to_step(sol)
        assert isinstance(result.step, np.ndarray)
        assert list(result.step) == [1, 2]


# ============== _EXPAND_STEP_TO_SOLUTION ==============

class TestExpandStepToSolution:
    """Tests for _expand_step_to_solution()."""

    def test_returns_solution_object(self):
        step = solutionstep()
        step.time = np.array([1.0, 2.0])
        step.Thickness = np.array([[1.0, 2.0], [3.0, 4.0]])
        result = _expand_step_to_solution(step)
        assert isinstance(result, solution)

    def test_number_of_steps_matches_time_length(self):
        step = solutionstep()
        step.time = np.array([1.0, 2.0, 3.0])
        step.Thickness = np.array([[1.0], [2.0], [3.0]])
        result = _expand_step_to_solution(step)
        assert len(result.steps) == 3

    def test_time_distributed_to_steps(self):
        step = solutionstep()
        step.time = np.array([0.5, 1.5])
        step.Vel = np.array([[10.0, 20.0], [30.0, 40.0]])
        result = _expand_step_to_solution(step)
        assert result.steps[0].time == pytest.approx(0.5)
        assert result.steps[1].time == pytest.approx(1.5)

    def test_array_split_along_time(self):
        step = solutionstep()
        step.time = np.array([1.0, 2.0])
        step.Vel = np.array([[1.0, 2.0], [3.0, 4.0]])
        result = _expand_step_to_solution(step)
        np.testing.assert_allclose(result.steps[0].Vel, [1.0, 2.0])
        np.testing.assert_allclose(result.steps[1].Vel, [3.0, 4.0])

    def test_single_scalar_time(self):
        step = solutionstep()
        step.time = 1.0  # scalar, not array
        result = _expand_step_to_solution(step)
        assert len(result.steps) == 1

    def test_no_time_attribute_falls_back(self):
        step = solutionstep()
        step.Thickness = np.array([[1.0, 2.0], [3.0, 4.0]])
        result = _expand_step_to_solution(step)
        # Should use first array dim (2) as nt
        assert len(result.steps) == 2


# ============== ROUND-TRIP ==============

class TestRoundTrip:
    """Test collapse -> expand round trip."""

    def test_round_trip_preserves_time(self):
        steps = [
            _make_step(1.0, 1, Thickness=np.array([10.0, 20.0])),
            _make_step(2.0, 2, Thickness=np.array([11.0, 21.0])),
        ]
        sol = _make_solution(steps)
        collapsed = _collapse_solution_to_step(sol)
        expanded = _expand_step_to_solution(collapsed)
        assert len(expanded.steps) == 2
        assert expanded.steps[0].time == pytest.approx(1.0)
        assert expanded.steps[1].time == pytest.approx(2.0)

    def test_round_trip_preserves_thickness(self):
        steps = [
            _make_step(1.0, 1, Thickness=np.array([10.0, 20.0])),
            _make_step(2.0, 2, Thickness=np.array([11.0, 21.0])),
        ]
        sol = _make_solution(steps)
        collapsed = _collapse_solution_to_step(sol)
        expanded = _expand_step_to_solution(collapsed)
        np.testing.assert_allclose(expanded.steps[0].Thickness, [10.0, 20.0])
        np.testing.assert_allclose(expanded.steps[1].Thickness, [11.0, 21.0])
