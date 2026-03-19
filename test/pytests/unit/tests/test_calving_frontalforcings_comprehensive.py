"""
Comprehensive unit tests for pyissm.model.classes.calving and frontalforcings modules.

Tests cover all calving and frontalforcings classes.
"""

import pytest

try:
    from pyissm.model.classes import calving, frontalforcings
    CLASSES_AVAILABLE = True
except ImportError:
    CLASSES_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not CLASSES_AVAILABLE,
    reason="Calving/frontalforcings classes not available"
)


# ============== CALVING TESTS ==============

class TestCalvingDefaultComprehensive:
    """Comprehensive tests for default calving class."""

    def test_init(self):
        c = calving.default()
        assert c is not None

    def test_has_calvingrate(self):
        c = calving.default()
        assert hasattr(c, 'calvingrate')

    def test_repr(self):
        c = calving.default()
        r = repr(c)
        assert isinstance(r, str)

    def test_str(self):
        c = calving.default()
        assert isinstance(str(c), str)


class TestCalvingVonmisesComprehensive:
    """Comprehensive tests for vonmises calving class."""

    def test_init(self):
        c = calving.vonmises()
        assert c is not None

    def test_has_stress_threshold_groundedice(self):
        c = calving.vonmises()
        assert hasattr(c, 'stress_threshold_groundedice')

    def test_has_stress_threshold_floatingice(self):
        c = calving.vonmises()
        assert hasattr(c, 'stress_threshold_floatingice')

    def test_has_min_thickness(self):
        c = calving.vonmises()
        assert hasattr(c, 'min_thickness')

    def test_repr(self):
        c = calving.vonmises()
        r = repr(c)
        assert isinstance(r, str)

    def test_str(self):
        c = calving.vonmises()
        assert isinstance(str(c), str)


class TestCalvingLevermannComprehensive:
    """Comprehensive tests for levermann calving class."""

    def test_init(self):
        c = calving.levermann()
        assert c is not None

    def test_has_coeff(self):
        c = calving.levermann()
        assert hasattr(c, 'coeff')

    def test_repr(self):
        c = calving.levermann()
        r = repr(c)
        assert isinstance(r, str)

    def test_str(self):
        c = calving.levermann()
        assert isinstance(str(c), str)


class TestCalvingMinthicknessComprehensive:
    """Comprehensive tests for minthickness calving class."""

    def test_init(self):
        c = calving.minthickness()
        assert c is not None

    def test_has_min_thickness(self):
        c = calving.minthickness()
        assert hasattr(c, 'min_thickness')

    def test_repr(self):
        c = calving.minthickness()
        r = repr(c)
        assert isinstance(r, str)

    def test_str(self):
        c = calving.minthickness()
        assert isinstance(str(c), str)


class TestCalvingCrevassedepthComprehensive:
    """Comprehensive tests for crevassedepth calving class."""

    def test_init(self):
        c = calving.crevassedepth()
        assert c is not None

    def test_has_crevasse_opening_stress(self):
        c = calving.crevassedepth()
        assert hasattr(c, 'crevasse_opening_stress')

    def test_has_water_height(self):
        c = calving.crevassedepth()
        assert hasattr(c, 'water_height')

    def test_repr(self):
        c = calving.crevassedepth()
        r = repr(c)
        assert isinstance(r, str)

    def test_str(self):
        c = calving.crevassedepth()
        assert isinstance(str(c), str)


class TestCalvingParameterizationComprehensive:
    """Comprehensive tests for parameterization calving class."""

    def test_init(self):
        c = calving.parameterization()
        assert c is not None

    def test_has_min_thickness(self):
        c = calving.parameterization()
        assert hasattr(c, 'min_thickness')

    def test_has_use_param(self):
        c = calving.parameterization()
        assert hasattr(c, 'use_param')

    def test_has_theta(self):
        c = calving.parameterization()
        assert hasattr(c, 'theta')

    def test_has_alpha(self):
        c = calving.parameterization()
        assert hasattr(c, 'alpha')

    def test_repr(self):
        c = calving.parameterization()
        r = repr(c)
        assert isinstance(r, str)

    def test_str(self):
        c = calving.parameterization()
        assert isinstance(str(c), str)


# ============== FRONTALFORCINGS TESTS ==============

class TestFrontalforcingsDefaultComprehensive:
    """Comprehensive tests for default frontalforcings class."""

    def test_init(self):
        f = frontalforcings.default()
        assert f is not None

    def test_has_meltingrate(self):
        f = frontalforcings.default()
        assert hasattr(f, 'meltingrate')

    def test_repr(self):
        f = frontalforcings.default()
        r = repr(f)
        assert isinstance(r, str)

    def test_str(self):
        f = frontalforcings.default()
        assert isinstance(str(f), str)


class TestFrontalforcingsRignotComprehensive:
    """Comprehensive tests for rignot frontalforcings class."""

    def test_init(self):
        f = frontalforcings.rignot()
        assert f is not None

    def test_has_thermalforcing(self):
        f = frontalforcings.rignot()
        assert hasattr(f, 'thermalforcing')

    def test_has_basin_id(self):
        f = frontalforcings.rignot()
        assert hasattr(f, 'basin_id')

    def test_has_num_basins(self):
        f = frontalforcings.rignot()
        assert hasattr(f, 'num_basins')
        assert f.num_basins == 0

    def test_has_subglacial_discharge(self):
        f = frontalforcings.rignot()
        assert hasattr(f, 'subglacial_discharge')

    def test_repr(self):
        f = frontalforcings.rignot()
        r = repr(f)
        assert isinstance(r, str)

    def test_str(self):
        f = frontalforcings.rignot()
        assert isinstance(str(f), str)


class TestFrontalforcingsRignotarmaComprehensive:
    """Comprehensive tests for rignotarma frontalforcings class."""

    def test_init(self):
        f = frontalforcings.rignotarma()
        assert f is not None

    def test_has_num_basins(self):
        f = frontalforcings.rignotarma()
        assert hasattr(f, 'num_basins')
        assert f.num_basins == 0

    def test_has_basin_id(self):
        f = frontalforcings.rignotarma()
        assert hasattr(f, 'basin_id')

    def test_has_ar_order(self):
        f = frontalforcings.rignotarma()
        assert hasattr(f, 'ar_order')

    def test_has_ma_order(self):
        f = frontalforcings.rignotarma()
        assert hasattr(f, 'ma_order')

    def test_has_num_breaks(self):
        f = frontalforcings.rignotarma()
        assert hasattr(f, 'num_breaks')
        assert f.num_breaks == 0

    def test_has_num_params(self):
        f = frontalforcings.rignotarma()
        assert hasattr(f, 'num_params')
        assert f.num_params == 0

    def test_repr(self):
        f = frontalforcings.rignotarma()
        r = repr(f)
        assert isinstance(r, str)

    def test_str(self):
        f = frontalforcings.rignotarma()
        assert isinstance(str(f), str)
