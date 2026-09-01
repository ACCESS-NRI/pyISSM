"""
Unit tests for pyissm.data.ocean_mask.gmtmask.

Property-based only (no old-vs-new comparison against the upstream MATLAB
gmtmask): the algorithm here is intentionally different (Natural Earth
polygons + vectorised point-in-polygon vs. a GMT CLI shell-out), so there
is no meaningful point-by-point "old vs. new" comparison to make. See
`goals.md` / `Implementation plan.md` for the rationale.
"""
from importlib import import_module

import numpy as np
import pytest

try:
    import cartopy
    CARTOPY_AVAILABLE = True
except ImportError:
    CARTOPY_AVAILABLE = False

pytestmark = pytest.mark.skipif(not CARTOPY_AVAILABLE, reason="cartopy not available")

if CARTOPY_AVAILABLE:
    ocean_mask_module = import_module('pyissm.data.ocean_mask')

RESOLUTION = '110m'


class TestGmtmaskKnownPoints:
    """Known ocean/land points, per the implementation plan's validation set."""

    def test_mid_pacific_is_ocean(self):
        mask = ocean_mask_module.gmtmask([0.], [-140.], resolution=RESOLUTION)
        assert mask[0] == 1

    def test_central_australia_is_land(self):
        mask = ocean_mask_module.gmtmask([-25.], [133.], resolution=RESOLUTION)
        assert mask[0] == 0

    def test_sahara_is_land(self):
        mask = ocean_mask_module.gmtmask([23.], [13.], resolution=RESOLUTION)
        assert mask[0] == 0


class TestGmtmaskShapeAndType:
    """Vectorised behaviour over multiple points at once."""

    def test_vectorised_call_shape_and_values(self):
        lat = np.array([0., -25., 23.])
        long = np.array([-140., 133., 13.])
        mask = ocean_mask_module.gmtmask(lat, long, resolution=RESOLUTION)
        assert mask.shape == (3,)
        assert np.array_equal(mask, [1, 0, 0])

    def test_returns_only_zeros_and_ones(self):
        lat = np.random.uniform(-89., 89., 200)
        long = np.random.uniform(-179., 179., 200)
        mask = ocean_mask_module.gmtmask(lat, long, resolution=RESOLUTION)
        assert set(np.unique(mask)).issubset({0, 1})


class TestGmtmaskValidation:

    def test_invalid_resolution_raises(self):
        with pytest.raises(ValueError):
            ocean_mask_module.gmtmask([0.], [0.], resolution='25m')
