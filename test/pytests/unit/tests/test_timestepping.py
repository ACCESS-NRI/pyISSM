"""
Unit tests for pyissm.model.classes.timestepping module.

Tests cover timestepping classes: default and adaptive.
"""

import pytest

try:
    from pyissm.model.classes import timestepping
    TIMESTEPPING_AVAILABLE = True
except ImportError:
    TIMESTEPPING_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not TIMESTEPPING_AVAILABLE,
    reason="timestepping classes not available"
)


class TestTimesteppingDefault:
    """Tests for default timestepping class."""

    def test_init_defaults(self):
        """Test default initialization."""
        t = timestepping.default()
        assert t is not None
        assert hasattr(t, 'start_time')
        assert hasattr(t, 'final_time')
        assert hasattr(t, 'time_step')

    def test_repr(self):
        """Test string representation."""
        t = timestepping.default()
        s = repr(t)
        assert isinstance(s, str)
        assert 'timestepping' in s.lower()

    def test_str(self):
        """Test short string."""
        t = timestepping.default()
        assert isinstance(str(t), str)

    def test_default_values(self):
        """Test default values are set correctly."""
        t = timestepping.default()
        assert t.start_time == 0
        # final_time may vary, just check it's a number
        assert isinstance(t.final_time, (int, float))
        assert t.time_step > 0


class TestTimesteppingAdaptive:
    """Tests for adaptive timestepping class."""

    def test_init_defaults(self):
        """Test default initialization."""
        t = timestepping.adaptive()
        assert t is not None
        assert hasattr(t, 'start_time')
        assert hasattr(t, 'final_time')
        assert hasattr(t, 'time_step_min')
        assert hasattr(t, 'time_step_max')
        assert hasattr(t, 'cfl_coefficient')

    def test_repr(self):
        """Test string representation."""
        t = timestepping.adaptive()
        s = repr(t)
        assert isinstance(s, str)
        assert 'adaptive' in s.lower() or 'timestepping' in s.lower()

    def test_str(self):
        """Test short string."""
        t = timestepping.adaptive()
        assert isinstance(str(t), str)

    def test_default_values(self):
        """Test default values are set correctly."""
        t = timestepping.adaptive()
        assert t.start_time == 0.0
        # Check CFL coefficient is reasonable
        assert 0.0 < t.cfl_coefficient <= 1.0
        # Check time step bounds
        assert t.time_step_min > 0
        assert t.time_step_max > t.time_step_min


class TestTimesteppingPhysics:
    """Tests for timestepping physical constraints."""

    def test_default_start_before_end(self):
        """Test default start time is before final time."""
        t = timestepping.default()
        assert t.start_time < t.final_time

    def test_adaptive_start_before_end(self):
        """Test adaptive start time is before final time."""
        t = timestepping.adaptive()
        assert t.start_time < t.final_time

    def test_default_step_positive(self):
        """Test default time step is positive."""
        t = timestepping.default()
        assert t.time_step > 0

    def test_adaptive_step_bounds_positive(self):
        """Test adaptive time step bounds are positive."""
        t = timestepping.adaptive()
        assert t.time_step_min > 0
        assert t.time_step_max > 0


class TestTimesteppingInheritance:
    """Tests for timestepping class inheritance behavior."""

    def test_default_from_other(self):
        """Test default initialization from another object."""
        t1 = timestepping.default()
        t1.start_time = 10.0
        t2 = timestepping.default(t1)
        assert t2.start_time == t1.start_time

    def test_adaptive_from_other(self):
        """Test adaptive initialization from another object."""
        t1 = timestepping.adaptive()
        t1.cfl_coefficient = 0.8
        t2 = timestepping.adaptive(t1)
        assert t2.cfl_coefficient == t1.cfl_coefficient
