"""
Unit tests for pyissm.model.classes.lovenumbers module.

Tests cover the packaged PREM Love-number lookup (`get_love_numbers`) and
the extended `lovenumbers` class `maxdeg`/`referenceframe` population.
"""

from importlib import import_module

import numpy as np
import pytest

try:
    # `pyissm.model.classes/__init__.py` re-exports the `lovenumbers` class
    # under the same name as this submodule, which shadows a plain
    # `import pyissm.model.classes.lovenumbers as x` (attribute lookup on
    # the package resolves to the class). Go through sys.modules instead.
    lovenumbers_module = import_module('pyissm.model.classes.lovenumbers')
    LOVENUMBERS_AVAILABLE = True
except ImportError:
    LOVENUMBERS_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not LOVENUMBERS_AVAILABLE,
    reason="lovenumbers module not available"
)


class TestGetLoveNumbers:
    """Tests for the get_love_numbers() lookup function."""

    def test_exact_match_known_prem_value(self):
        """Exact-match check against a known PREM value (the one function validated this way)."""
        series = lovenumbers_module.get_love_numbers('loadingverticaldisplacement', maxdeg=10)
        assert len(series) == 11
        assert series[0] == 0
        assert series[1] == pytest.approx(-1.28740059, abs=1e-6)

    def test_cf_reference_frame_overrides_degree_1(self):
        """CF reference frame applies the Blewitt (2003) degree-1 override."""
        cm_series = lovenumbers_module.get_love_numbers('loadingverticaldisplacement', referenceframe='CM', maxdeg=10)
        cf_series = lovenumbers_module.get_love_numbers('loadingverticaldisplacement', referenceframe='CF', maxdeg=10)
        assert cf_series[1] == pytest.approx(-0.269)
        assert cm_series[1] != pytest.approx(-0.269)
        # Higher degrees are unaffected by the degree-1 override
        assert cf_series[2] == cm_series[2]

    def test_cf_override_only_applies_to_defined_types(self):
        """Tidal love-number types have no CF degree-1 override defined."""
        cm_series = lovenumbers_module.get_love_numbers('tidalverticaldisplacement', referenceframe='CM', maxdeg=10)
        cf_series = lovenumbers_module.get_love_numbers('tidalverticaldisplacement', referenceframe='CF', maxdeg=10)
        assert cf_series[1] == cm_series[1]

    def test_invalid_love_type_raises(self):
        """An unrecognized love_type raises ValueError."""
        with pytest.raises(ValueError):
            lovenumbers_module.get_love_numbers('not_a_real_type', maxdeg=10)

    def test_invalid_referenceframe_raises(self):
        """A referenceframe other than 'CM'/'CF' raises ValueError."""
        with pytest.raises(ValueError):
            lovenumbers_module.get_love_numbers('loadingverticaldisplacement', referenceframe='XX', maxdeg=10)

    def test_maxdeg_exceeds_table_raises(self):
        """maxdeg above the tabulated 10000 degrees raises ValueError."""
        with pytest.raises(ValueError):
            lovenumbers_module.get_love_numbers('loadingverticaldisplacement', maxdeg=10001)


class TestLovenumbersClass:
    """Tests for the lovenumbers class, including maxdeg-based population."""

    def test_default_init_yields_all_nan(self):
        """Regression: lovenumbers() with no maxdeg still yields the original all-NaN defaults."""
        ln = lovenumbers_module.lovenumbers()
        assert np.isnan(ln.h)
        assert np.isnan(ln.k)
        assert np.isnan(ln.l)
        assert np.isnan(ln.th)
        assert np.isnan(ln.tk)
        assert np.isnan(ln.tl)
        assert np.isnan(ln.pmtf_colinear)
        assert np.isnan(ln.pmtf_ortho)
        assert ln.istime == 1

    def test_maxdeg_populates_all_fields(self):
        """lovenumbers(maxdeg=100) populates h/k/l/th/tk/tl as (101, 1) arrays with no NaN."""
        ln = lovenumbers_module.lovenumbers(maxdeg=100)
        for field in ('h', 'k', 'l', 'th', 'tk', 'tl'):
            array = getattr(ln, field)
            assert array.shape == (101, 1)
            assert not np.any(np.isnan(array))
        assert np.all(np.isfinite(ln.pmtf_colinear))
        assert np.all(np.isfinite(ln.pmtf_ortho))
        assert ln.timefreq.tolist() == [0]
        assert ln.istime == 1

    def test_maxdeg_matches_get_love_numbers(self):
        """The populated h array matches get_love_numbers() directly (exact-match target)."""
        ln = lovenumbers_module.lovenumbers(maxdeg=10)
        expected = lovenumbers_module.get_love_numbers('loadingverticaldisplacement', 'CM', 10)
        np.testing.assert_array_equal(ln.h.reshape(-1), expected)

    def test_maxdeg_below_degree_2_leaves_pmtf_default(self):
        """pmtf_colinear/pmtf_ortho need degree 2 and stay at [[0.0]] below it."""
        ln = lovenumbers_module.lovenumbers(maxdeg=1)
        assert ln.pmtf_colinear.tolist() == [[0.0]]
        assert ln.pmtf_ortho.tolist() == [[0.0]]

    def test_maxdeg_takes_precedence_over_other(self):
        """An explicit maxdeg= wins over a NaN-filled `other` instance."""
        nan_instance = lovenumbers_module.lovenumbers()
        ln = lovenumbers_module.lovenumbers(other=nan_instance, maxdeg=10)
        assert not np.any(np.isnan(ln.h))

    def test_repr_no_longer_reports_not_implemented(self):
        """The stale 'NOT YET IMPLEMENTED' banner has been removed from __repr__."""
        ln = lovenumbers_module.lovenumbers(maxdeg=10)
        assert 'NOT YET IMPLEMENTED' not in repr(ln)
