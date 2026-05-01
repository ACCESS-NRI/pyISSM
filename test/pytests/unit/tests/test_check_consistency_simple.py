"""
Unit tests for check_consistency() in simple model classes (timestepping.default).
Uses SimpleNamespace mock for md.
"""

import pytest
import numpy as np
from types import SimpleNamespace

try:
    from pyissm.model.classes.timestepping import default as timestepping_default
    TIMESTEPPING_AVAILABLE = True
except ImportError:
    TIMESTEPPING_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not TIMESTEPPING_AVAILABLE,
    reason="pyissm.model.classes.timestepping not available"
)


def _make_md():
    """Create a minimal mock md for check_consistency."""
    md = SimpleNamespace()
    md.mesh = SimpleNamespace()
    md.mesh.numberofvertices = 10
    md.mesh.numberofelements = 8
    md.private = SimpleNamespace()
    md.private.isconsistent = True
    messages = []
    def check_message(s):
        md.private.isconsistent = False
        messages.append(s)
    md.check_message = check_message
    md._messages = messages
    return md


def _attach_timestepping(md, ts):
    """Attach a timestepping object to md so _check_field can resolve 'timestepping.X'."""
    md.timestepping = ts


class TestTimesteppingDefaultCheckConsistency:
    """Tests for timestepping.default.check_consistency()."""

    def test_valid_defaults_consistent(self):
        ts = timestepping_default()
        md = _make_md()
        _attach_timestepping(md, ts)
        ts.check_consistency(md, solution='', analyses=[])
        assert md.private.isconsistent is True

    def test_returns_md(self):
        ts = timestepping_default()
        md = _make_md()
        _attach_timestepping(md, ts)
        result = ts.check_consistency(md, solution='', analyses=[])
        assert result is md

    def test_negative_time_range_fails(self):
        ts = timestepping_default()
        ts.start_time = 10
        ts.final_time = 5  # final < start
        md = _make_md()
        _attach_timestepping(md, ts)
        ts.check_consistency(md, solution='', analyses=[])
        assert md.private.isconsistent is False

    def test_equal_start_final_time_passes(self):
        ts = timestepping_default()
        ts.start_time = 5
        ts.final_time = 5  # equal (range = 0, not negative)
        md = _make_md()
        _attach_timestepping(md, ts)
        ts.check_consistency(md, solution='', analyses=[])
        assert md.private.isconsistent is True

    def test_invalid_interp_forcing_fails(self):
        ts = timestepping_default()
        ts.interp_forcing = 5  # only 0 or 1 allowed
        md = _make_md()
        _attach_timestepping(md, ts)
        ts.check_consistency(md, solution='', analyses=[])
        assert md.private.isconsistent is False

    def test_invalid_average_forcing_fails(self):
        ts = timestepping_default()
        ts.average_forcing = 2  # only 0 or 1 allowed
        md = _make_md()
        _attach_timestepping(md, ts)
        ts.check_consistency(md, solution='', analyses=[])
        assert md.private.isconsistent is False

    def test_invalid_cycle_forcing_fails(self):
        ts = timestepping_default()
        ts.cycle_forcing = -1  # only 0 or 1 allowed
        md = _make_md()
        _attach_timestepping(md, ts)
        ts.check_consistency(md, solution='', analyses=[])
        assert md.private.isconsistent is False

    def test_nan_start_time_fails(self):
        ts = timestepping_default()
        ts.start_time = np.nan
        md = _make_md()
        _attach_timestepping(md, ts)
        ts.check_consistency(md, solution='', analyses=[])
        assert md.private.isconsistent is False

    def test_nan_final_time_fails(self):
        ts = timestepping_default()
        ts.final_time = np.nan
        md = _make_md()
        _attach_timestepping(md, ts)
        ts.check_consistency(md, solution='', analyses=[])
        assert md.private.isconsistent is False

    def test_nan_time_step_fails(self):
        ts = timestepping_default()
        ts.time_step = np.nan
        md = _make_md()
        _attach_timestepping(md, ts)
        ts.check_consistency(md, solution='', analyses=[])
        assert md.private.isconsistent is False

    def test_negative_time_step_fails(self):
        ts = timestepping_default()
        ts.time_step = -1.0
        md = _make_md()
        _attach_timestepping(md, ts)
        ts.check_consistency(md, solution='', analyses=[])
        assert md.private.isconsistent is False

    def test_zero_time_step_passes_non_transient(self):
        ts = timestepping_default()
        ts.time_step = 0.0  # ge=0 is ok for non-transient solution
        md = _make_md()
        _attach_timestepping(md, ts)
        ts.check_consistency(md, solution='StressbalanceSolution', analyses=[])
        assert md.private.isconsistent is True

    def test_zero_time_step_fails_transient(self):
        ts = timestepping_default()
        ts.time_step = 0.0  # gt=0 required for TransientSolution
        md = _make_md()
        _attach_timestepping(md, ts)
        ts.check_consistency(md, solution='TransientSolution', analyses=[])
        assert md.private.isconsistent is False

    def test_positive_time_step_passes_transient(self):
        ts = timestepping_default()
        ts.time_step = 0.5
        md = _make_md()
        _attach_timestepping(md, ts)
        ts.check_consistency(md, solution='TransientSolution', analyses=[])
        assert md.private.isconsistent is True

    def test_large_valid_time_range(self):
        ts = timestepping_default()
        ts.start_time = 0
        ts.final_time = 1000
        ts.time_step = 1.0
        md = _make_md()
        _attach_timestepping(md, ts)
        ts.check_consistency(md, solution='TransientSolution', analyses=[])
        assert md.private.isconsistent is True
