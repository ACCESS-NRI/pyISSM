"""
Comprehensive unit tests for pyissm.model.classes.hydrology module.

Tests cover all hydrology classes.
"""

import pytest

try:
    from pyissm.model.classes import hydrology
    HYDRO_AVAILABLE = True
except ImportError:
    HYDRO_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not HYDRO_AVAILABLE,
    reason="Hydrology classes not available"
)


class TestHydrologyShreveComprehensive:
    """Comprehensive tests for shreve hydrology class."""

    def test_init(self):
        h = hydrology.shreve()
        assert h is not None

    def test_has_spcwatercolumn(self):
        h = hydrology.shreve()
        assert hasattr(h, 'spcwatercolumn')

    def test_has_stabilization(self):
        h = hydrology.shreve()
        assert hasattr(h, 'stabilization')
        assert h.stabilization == 1

    def test_repr(self):
        h = hydrology.shreve()
        r = repr(h)
        assert isinstance(r, str)

    def test_str(self):
        h = hydrology.shreve()
        assert isinstance(str(h), str)


class TestHydrologyDcComprehensive:
    """Comprehensive tests for dc hydrology class."""

    def test_init(self):
        h = hydrology.dc()
        assert h is not None

    def test_has_water_compressibility(self):
        h = hydrology.dc()
        assert hasattr(h, 'water_compressibility')

    def test_has_isefficientlayer(self):
        h = hydrology.dc()
        assert hasattr(h, 'isefficientlayer')
        assert h.isefficientlayer == 1

    def test_has_penalty_factor(self):
        h = hydrology.dc()
        assert hasattr(h, 'penalty_factor')
        assert h.penalty_factor == 3

    def test_has_penalty_lock(self):
        h = hydrology.dc()
        assert hasattr(h, 'penalty_lock')
        assert h.penalty_lock == 0

    def test_has_rel_tol(self):
        h = hydrology.dc()
        assert hasattr(h, 'rel_tol')

    def test_has_max_iter(self):
        h = hydrology.dc()
        assert hasattr(h, 'max_iter')
        assert h.max_iter == 100

    def test_has_sediment_thickness(self):
        h = hydrology.dc()
        assert hasattr(h, 'sediment_thickness')

    def test_has_sediment_transmitivity(self):
        h = hydrology.dc()
        assert hasattr(h, 'sediment_transmitivity')

    def test_has_epl_conductivity(self):
        h = hydrology.dc()
        assert hasattr(h, 'epl_conductivity')

    def test_repr(self):
        h = hydrology.dc()
        r = repr(h)
        assert isinstance(r, str)

    def test_str(self):
        h = hydrology.dc()
        assert isinstance(str(h), str)


class TestHydrologyGladsComprehensive:
    """Comprehensive tests for glads hydrology class."""

    def test_init(self):
        h = hydrology.glads()
        assert h is not None

    def test_has_pressure_melt_coefficient(self):
        h = hydrology.glads()
        assert hasattr(h, 'pressure_melt_coefficient')

    def test_has_sheet_conductivity(self):
        h = hydrology.glads()
        assert hasattr(h, 'sheet_conductivity')

    def test_has_cavity_spacing(self):
        h = hydrology.glads()
        assert hasattr(h, 'cavity_spacing')
        assert h.cavity_spacing == 2.0

    def test_has_bump_height(self):
        h = hydrology.glads()
        assert hasattr(h, 'bump_height')

    def test_has_ischannels(self):
        h = hydrology.glads()
        assert hasattr(h, 'ischannels')
        assert h.ischannels == False

    def test_has_channel_conductivity(self):
        h = hydrology.glads()
        assert hasattr(h, 'channel_conductivity')

    def test_has_englacial_void_ratio(self):
        h = hydrology.glads()
        assert hasattr(h, 'englacial_void_ratio')

    def test_repr(self):
        h = hydrology.glads()
        r = repr(h)
        assert isinstance(r, str)

    def test_str(self):
        h = hydrology.glads()
        assert isinstance(str(h), str)


class TestHydrologyShaktiComprehensive:
    """Comprehensive tests for shakti hydrology class."""

    def test_init(self):
        h = hydrology.shakti()
        assert h is not None

    def test_has_spchead(self):
        h = hydrology.shakti()
        assert hasattr(h, 'spchead')

    def test_has_relaxation(self):
        h = hydrology.shakti()
        assert hasattr(h, 'relaxation')

    def test_has_storage(self):
        h = hydrology.shakti()
        assert hasattr(h, 'storage')

    def test_repr(self):
        h = hydrology.shakti()
        r = repr(h)
        assert isinstance(r, str)

    def test_str(self):
        h = hydrology.shakti()
        assert isinstance(str(h), str)


class TestHydrologyPismComprehensive:
    """Comprehensive tests for pism hydrology class."""

    def test_init(self):
        h = hydrology.pism()
        assert h is not None

    def test_has_drainage_rate(self):
        h = hydrology.pism()
        assert hasattr(h, 'drainage_rate')

    def test_has_watercolumn_max(self):
        h = hydrology.pism()
        assert hasattr(h, 'watercolumn_max')

    def test_has_requested_outputs(self):
        h = hydrology.pism()
        assert hasattr(h, 'requested_outputs')

    def test_repr(self):
        h = hydrology.pism()
        r = repr(h)
        assert isinstance(r, str)

    def test_str(self):
        h = hydrology.pism()
        assert isinstance(str(h), str)


class TestHydrologyArmapwComprehensive:
    """Comprehensive tests for armapw hydrology class."""

    def test_init(self):
        h = hydrology.armapw()
        assert h is not None

    def test_has_num_basins(self):
        h = hydrology.armapw()
        assert hasattr(h, 'num_basins')
        assert h.num_basins == 0

    def test_has_basin_id(self):
        h = hydrology.armapw()
        assert hasattr(h, 'basin_id')

    def test_has_ar_order(self):
        h = hydrology.armapw()
        assert hasattr(h, 'ar_order')

    def test_has_ma_order(self):
        h = hydrology.armapw()
        assert hasattr(h, 'ma_order')

    def test_repr(self):
        h = hydrology.armapw()
        r = repr(h)
        assert isinstance(r, str)

    def test_str(self):
        h = hydrology.armapw()
        assert isinstance(str(h), str)


class TestHydrologyTwsComprehensive:
    """Comprehensive tests for tws hydrology class."""

    def test_init(self):
        h = hydrology.tws()
        assert h is not None

    def test_has_spcwatercolumn(self):
        h = hydrology.tws()
        assert hasattr(h, 'spcwatercolumn')

    def test_repr(self):
        h = hydrology.tws()
        r = repr(h)
        assert isinstance(r, str)

    def test_str(self):
        h = hydrology.tws()
        assert isinstance(str(h), str)
