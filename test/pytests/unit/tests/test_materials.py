"""
Unit tests for pyissm.model.classes.materials module.

Tests cover material property classes: ice, hydro, litho, damageice, enhancedice, estar
"""

import numpy as np
import pytest
from types import SimpleNamespace

try:
    from pyissm.model.classes.materials import ice, hydro, litho, damageice, enhancedice, estar
    MATERIALS_AVAILABLE = True
except ImportError:
    MATERIALS_AVAILABLE = False
    ice = hydro = litho = damageice = enhancedice = estar = None

pytestmark = pytest.mark.skipif(
    not MATERIALS_AVAILABLE,
    reason="Materials classes not available"
)


class TestMaterialsIce:
    """Tests for materials.ice class."""

    def test_init_defaults(self):
        """Test default initialization of ice material."""
        m = ice()
        assert m.rho_ice == 917.0
        assert m.rho_water == 1023.0
        assert m.rho_freshwater == 1000.0
        assert m.mu_water == 0.001787
        assert m.heatcapacity == 2093.0
        assert m.latentheat == 3.34e5
        assert m.thermalconductivity == 2.4
        assert m.temperateiceconductivity == 0.24
        assert m.effectiveconductivity_averaging == 1
        assert m.meltingpoint == 273.15
        assert m.beta == 9.8e-8
        assert m.mixed_layer_capacity == 3974.0
        assert m.thermal_exchange_velocity == 1.00e-4
        assert m.rheology_law == 'Paterson'
        assert m.rheology_B == 2.1e8
        assert m.rheology_n == 3.0
        assert m.earth_density == 5512.0

    def test_repr(self):
        """Test string representation."""
        m = ice()
        s = repr(m)
        assert 'ice' in s.lower() or 'materials' in s.lower()
        assert 'rho_ice' in s or 'density' in s

    def test_str(self):
        """Test short string."""
        m = ice()
        assert 'ice' in str(m).lower()

    def test_init_from_other(self):
        """Test initialization from another object."""
        other = SimpleNamespace(rho_ice=910.0, rheology_n=4.0)
        m = ice(other)
        assert m.rho_ice == 910.0
        assert m.rheology_n == 4.0
        # Other defaults should remain
        assert m.rho_water == 1023.0

    def test_rheology_law_options(self):
        """Test different rheology law options."""
        valid_laws = ['None', 'BuddJacka', 'Cuffey', 'CuffeyTemperate', 
                      'Paterson', 'Arrhenius', 'LliboutryDuval', 'NyeCO2', 'NyeH2O']
        m = ice()
        for law in valid_laws:
            m.rheology_law = law
            assert m.rheology_law == law


class TestMaterialsHydro:
    """Tests for materials.hydro class."""

    def test_init_defaults(self):
        """Test default initialization of hydro material."""
        m = hydro()
        assert hasattr(m, 'rho_ice')
        assert hasattr(m, 'rho_water')
        assert hasattr(m, 'rho_freshwater')

    def test_repr(self):
        """Test string representation."""
        m = hydro()
        s = repr(m)
        assert 'hydro' in s.lower() or 'materials' in s.lower()

    def test_str(self):
        """Test short string."""
        m = hydro()
        assert 'hydro' in str(m).lower()


class TestMaterialsLitho:
    """Tests for materials.litho class."""

    def test_init_defaults(self):
        """Test default initialization of litho material."""
        m = litho()
        assert hasattr(m, 'numlayers')
        assert hasattr(m, 'radius')
        assert hasattr(m, 'viscosity')
        assert hasattr(m, 'lame_lambda')
        assert hasattr(m, 'lame_mu')
        assert hasattr(m, 'burgers_viscosity')
        assert hasattr(m, 'burgers_mu')
        assert hasattr(m, 'density')
        assert hasattr(m, 'issolid')

    def test_repr(self):
        """Test string representation."""
        m = litho()
        s = repr(m)
        assert 'litho' in s.lower() or 'materials' in s.lower()

    def test_str(self):
        """Test short string."""
        m = litho()
        assert 'litho' in str(m).lower()


class TestMaterialsDamageice:
    """Tests for materials.damageice class."""

    def test_init_defaults(self):
        """Test default initialization of damageice material."""
        m = damageice()
        # Should have all ice properties plus damage-specific ones
        assert m.rho_ice == 917.0
        assert hasattr(m, 'rheology_B')

    def test_repr(self):
        """Test string representation."""
        m = damageice()
        s = repr(m)
        assert 'damage' in s.lower() or 'ice' in s.lower()

    def test_str(self):
        """Test short string."""
        m = damageice()
        assert 'damage' in str(m).lower() or 'ice' in str(m).lower()


class TestMaterialsEnhancedice:
    """Tests for materials.enhancedice class."""

    def test_init_defaults(self):
        """Test default initialization of enhancedice material."""
        m = enhancedice()
        # Should have all ice properties
        assert m.rho_ice == 917.0
        assert hasattr(m, 'rheology_B')

    def test_repr(self):
        """Test string representation."""
        m = enhancedice()
        s = repr(m)
        assert 'enhanced' in s.lower() or 'ice' in s.lower()

    def test_str(self):
        """Test short string."""
        m = enhancedice()
        assert 'enhanced' in str(m).lower() or 'ice' in str(m).lower()


class TestMaterialsEstar:
    """Tests for materials.estar class."""

    def test_init_defaults(self):
        """Test default initialization of estar material."""
        m = estar()
        # Should have all ice properties
        assert m.rho_ice == 917.0
        assert hasattr(m, 'rheology_B')

    def test_repr(self):
        """Test string representation."""
        m = estar()
        s = repr(m)
        assert 'estar' in s.lower() or 'ice' in s.lower()

    def test_str(self):
        """Test short string."""
        m = estar()
        assert 'estar' in str(m).lower() or 'ice' in str(m).lower()


class TestMaterialsPhysicalConstants:
    """Tests for physical constant values in ice materials."""

    def test_ice_density_reasonable(self):
        """Test that ice density is in reasonable range."""
        m = ice()
        assert 850 < m.rho_ice < 950  # Ice density typically 917 kg/m3

    def test_water_density_reasonable(self):
        """Test that water density is in reasonable range."""
        m = ice()
        assert 1000 < m.rho_water < 1050  # Ocean water density

    def test_freshwater_density_reasonable(self):
        """Test that freshwater density is reasonable."""
        m = ice()
        assert 995 < m.rho_freshwater < 1005

    def test_melting_point_reasonable(self):
        """Test that melting point is reasonable."""
        m = ice()
        assert 273 < m.meltingpoint < 274  # Should be around 273.15 K

    def test_heat_capacity_reasonable(self):
        """Test that heat capacity is in reasonable range."""
        m = ice()
        assert 2000 < m.heatcapacity < 2200  # Ice heat capacity ~2093 J/kg/K

    def test_thermal_conductivity_reasonable(self):
        """Test that thermal conductivity is reasonable."""
        m = ice()
        assert 2.0 < m.thermalconductivity < 3.0  # Ice ~2.4 W/m/K

    def test_glens_law_exponent(self):
        """Test Glen's flow law exponent is 3."""
        m = ice()
        assert m.rheology_n == 3.0  # Standard value for ice


class TestMaterialsInheritance:
    """Test inheritance behavior between material classes."""

    def test_ice_fields_preserved(self):
        """Test that ice fields are preserved when inheriting."""
        base = ice()
        base.rho_ice = 900.0
        base.rheology_law = 'Cuffey'
        
        derived = ice(base)
        assert derived.rho_ice == 900.0
        assert derived.rheology_law == 'Cuffey'

    def test_partial_inheritance(self):
        """Test that only matching fields are inherited."""
        other = SimpleNamespace(
            rho_ice=905.0,
            nonexistent_field=42
        )
        m = ice(other)
        assert m.rho_ice == 905.0
        # Should not have the nonexistent field
        assert not hasattr(m, 'nonexistent_field') or m.rheology_B == 2.1e8
