"""
Unit tests for pyissm.tools.materials - ice rigidity functions.
"""

import pytest
import warnings
import numpy as np

try:
    from pyissm.tools.materials import paterson, cuffey, nye
    MATERIALS_AVAILABLE = True
except ImportError:
    MATERIALS_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not MATERIALS_AVAILABLE,
    reason="pyissm.tools.materials not available"
)


# ============== CUFFEY TESTS ==============

class TestCuffey:
    """Tests for the cuffey() ice rigidity function."""

    def test_scalar_input_returns_array(self):
        result = cuffey(np.array(263.15))
        assert isinstance(result, np.ndarray)

    def test_cold_ice_positive_rigidity(self):
        # Cold ice (~-30 C = 243.15 K) should have positive rigidity
        result = cuffey(np.array([243.15]))
        assert result[0] > 0

    def test_near_melting_positive_rigidity(self):
        # Near-melting ice (~-2 C = 271.15 K) should have positive rigidity
        result = cuffey(np.array([271.15]))
        assert result[0] > 0

    def test_array_input(self):
        temps = np.array([220.0, 240.0, 260.0, 270.0])
        result = cuffey(temps)
        assert result.shape == temps.shape

    def test_output_all_positive(self):
        temps = np.linspace(200.0, 272.0, 50)
        result = cuffey(temps)
        assert np.all(result > 0)

    def test_warmer_ice_lower_rigidity(self):
        # Warmer ice should be less rigid (lower B value)
        cold = cuffey(np.array([230.0]))[0]
        warm = cuffey(np.array([265.0]))[0]
        assert cold > warm

    def test_negative_temperature_raises(self):
        with pytest.raises(RuntimeError, match="Kelvin"):
            cuffey(np.array([-10.0]))

    def test_zero_kelvin_raises(self):
        with pytest.raises(RuntimeError):
            cuffey(np.array([0.0]))

    def test_floor_value_applied(self):
        # Very warm temperatures might compute negative rigidity before floor
        # Ensure floor of 1e6 is applied
        result = cuffey(np.array([272.0]))
        assert result[0] >= 1e6

    def test_returns_float_dtype(self):
        result = cuffey(np.array([250.0]))
        assert result.dtype == float

    def test_shape_preserved_2d(self):
        temps = np.array([[230.0, 240.0], [250.0, 260.0]])
        result = cuffey(temps)
        assert result.shape == temps.shape

    def test_piecewise_boundaries(self):
        # Test temperatures near piecewise boundaries
        boundaries = np.array([228.15, 233.15, 238.15, 243.15, 248.15, 253.15, 258.15, 263.15, 268.15])
        result = cuffey(boundaries)
        assert np.all(result > 0)
        assert np.all(np.isfinite(result))


# ============== PATERSON TESTS ==============

class TestPaterson:
    """Tests for the paterson() ice rigidity function (deprecated)."""

    def test_issues_deprecation_warning(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            paterson(np.array([250.0]))
            assert len(w) >= 1
            assert any(issubclass(warning.category, DeprecationWarning) for warning in w)

    def test_returns_array(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = paterson(np.array([250.0]))
        assert isinstance(result, np.ndarray)

    def test_positive_rigidity(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = paterson(np.array([240.0, 250.0, 260.0]))
        assert np.all(result > 0)

    def test_negative_temperature_raises(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            with pytest.raises(RuntimeError, match="Kelvin"):
                paterson(np.array([-5.0]))

    def test_array_input(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            temps = np.linspace(200.0, 272.0, 30)
            result = paterson(temps)
        assert result.shape == temps.shape
        assert np.all(result > 0)

    def test_floor_applied(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = paterson(np.array([272.0]))
        assert result[0] >= 1e6

    def test_piecewise_segments_covered(self):
        # Test temperatures in each piecewise segment
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            temps = np.array([220.0, 227.0, 232.0, 237.0, 242.0, 247.0, 252.0, 257.0, 262.0, 267.0, 271.0])
            result = paterson(temps)
        assert np.all(result > 0)
        assert np.all(np.isfinite(result))

    def test_returns_float_dtype(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = paterson(np.array([250.0]))
        assert result.dtype == float


# ============== NYE TESTS ==============

class TestNye:
    """Tests for the nye() ice rigidity function."""

    # ---- H2O ice tests ----

    def test_h2o_scalar_returns_array(self):
        result = nye(np.array([250.0]), ice_type=2)
        assert isinstance(result, np.ndarray)

    def test_h2o_cold_positive_rigidity(self):
        result = nye(np.array([200.0]), ice_type=2)
        assert result[0] > 0

    def test_h2o_array_input(self):
        temps = np.array([200.0, 220.0, 240.0, 260.0])
        result = nye(temps, ice_type=2)
        assert result.shape == temps.shape

    def test_h2o_all_positive(self):
        temps = np.linspace(200.0, 272.0, 20)
        result = nye(temps, ice_type=2)
        assert np.all(result > 0)

    def test_h2o_melting_warning(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            nye(np.array([280.0]), ice_type=2)
            assert any("melting" in str(warning.message).lower() for warning in w)

    # ---- CO2 ice tests ----

    def test_co2_returns_array(self):
        result = nye(np.array([150.0]), ice_type=1)
        assert isinstance(result, np.ndarray)

    def test_co2_cold_positive_rigidity(self):
        result = nye(np.array([150.0]), ice_type=1)
        assert result[0] > 0

    def test_co2_array_input(self):
        temps = np.array([100.0, 150.0, 180.0])
        result = nye(temps, ice_type=1)
        assert result.shape == temps.shape

    def test_co2_possible_melting_warning(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            nye(np.array([210.0]), ice_type=1)
            assert any("possible melting" in str(warning.message).lower() for warning in w)

    def test_co2_guaranteed_melting_warning(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            nye(np.array([225.0]), ice_type=1)
            assert any("guaranteed melting" in str(warning.message).lower() for warning in w)

    # ---- Error cases ----

    def test_invalid_ice_type_raises(self):
        with pytest.raises(ValueError, match="ice_type must be"):
            nye(np.array([250.0]), ice_type=3)

    def test_invalid_ice_type_zero_raises(self):
        with pytest.raises(ValueError):
            nye(np.array([250.0]), ice_type=0)

    def test_h2o_warmer_ice_lower_rigidity(self):
        # Warmer H2O ice should generally be less rigid
        cold = nye(np.array([200.0]), ice_type=2)[0]
        warm = nye(np.array([260.0]), ice_type=2)[0]
        assert cold > warm
