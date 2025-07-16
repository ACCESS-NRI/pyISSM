import numpy as np
from . import build_utils
from . import class_registry

## ------------------------------------------------------
## materials.ice
## ------------------------------------------------------
@class_registry.register_class
class ice(class_registry.manage_state):
    '''
    materials.ice Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.rho_ice = 917.
        self.rho_water = 1023.
        self.rho_freshwater = 1000.
        self.mu_water = 0.001787
        self.heatcapacity = 2093.
        self.latentheat = 3.34e5
        self.thermalconductivity = 2.4
        self.temperateiceconductivity = 0.24
        self.effectiveconductivity_averaging = 1
        self.meltingpoint = 273.15
        self.beta = 9.8e-8
        self.mixed_layer_capacity = 3974.
        self.thermal_exchange_velocity = 1.00e-4
        self.rheology_law = 'Paterson'
        self.rheology_B = 2.1 * 1e8
        self.rheology_n = 3.
        self.earth_density = 5512.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Materials (ice):\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'rho_ice', 'ice density [kg/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rho_water', 'ocean water density [kg/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rho_freshwater', 'fresh water density [kg/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'mu_water', 'water viscosity [N s/m^2]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'heatcapacity', 'heat capacity [J/kg/K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thermalconductivity', 'ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'temperateiceconductivity', 'temperate ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'meltingpoint', 'melting point of ice at 1atm in K'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'latentheat', 'latent heat of fusion [J/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'beta', 'rate of change of melting point with pressure [K/Pa]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'mixed_layer_capacity', 'mixed layer capacity [W/kg/K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thermal_exchange_velocity', 'thermal exchange velocity [m/s]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rheology_B', 'flow law parameter [Pa s^(1/n)]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rheology_n', 'Glen\'s flow law exponent'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rheology_law', 'law for the temperature dependance of the rheology: \'None\', \'BuddJacka\', \'Cuffey\', \'CuffeyTemperate\', \'Paterson\', \'Arrhenius\', \'LliboutryDuval\', \'NyeCO2\', or \'NyeH2O\''))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - materials.ice Class'
        return s

## ------------------------------------------------------
## materials.hydro
## ------------------------------------------------------
@class_registry.register_class
class hydro(class_registry.manage_state):
    '''
    materials.hydro Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.rho_ice = 917.
        self.rho_water = 1023.
        self.rho_freshwater = 1000.
        self.earth_density = 5512.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Materials (hydro):\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'rho_ice', 'ice density [kg/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rho_water', 'ocean water density [kg/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rho_freshwater', 'fresh water density [kg/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'earth_density', 'mantle density [kg/m^3]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - materials.hydro Class'
        return s

## ------------------------------------------------------
## materials.litho
## ------------------------------------------------------
@class_registry.register_class
class litho(class_registry.manage_state):
    '''
    materials.litho Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.numlayers = 2.
        self.radius = [1e3, 6278e3, 6378e3]
        self.viscosity = [1e21, 1e40]
        self.lame_mu = [1.45e11, 6.7e10]
        self.lame_lambda = self.lame_mu
        self.burgers_viscosity = [np.nan, np.nan]
        self.burgers_mu = [np.nan, np.nan]
        self.ebm_alpha = [np.nan, np.nan]
        self.ebm_delta = [np.nan, np.nan]
        self.ebm_taul = [np.nan, np.nan]
        self.ebm_tauh = [np.nan, np.nan]
        self.rheologymodel = [0, 0]
        self.density = [5.51e3, 5.50e3]
        self.issolid = [1, 1]
        self.earth_density = 5512.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Materials (litho):\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'numlayers', 'number of layers (default: 2)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'radius', 'array describing the radius for each interface (numlayers + 1) [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'viscosity', 'array describing each layer\'s viscosity (numlayers) [Pa.s]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'lame_lambda', 'array describing the lame lambda parameter (numlayers) [Pa]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'lame_mu', 'array describing the shear modulus for each layers (numlayers) [Pa]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'burgers_viscosity', 'array describing each layer\'s transient viscosity, only for Burgers rheologies  (numlayers) [Pa.s]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'burgers_mu', 'array describing each layer\'s transient shear modulus, only for Burgers rheologies  (numlayers) [Pa]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'ebm_alpha', 'array describing each layer\'s exponent parameter controlling the shape of shear modulus curve between taul and tauh, only for EBM rheology (numlayers)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'ebm_delta', 'array describing each layer\'s amplitude of the transient relaxation (ratio between elastic rigity to pre-maxwell relaxation rigity), only for EBM rheology (numlayers)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'ebm_taul', 'array describing each layer\'s starting period for transient relaxation, only for EBM rheology  (numlayers) [s]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'ebm_tauh', 'array describing each layer''s array describing each layer\'s end period for transient relaxation, only for Burgers rheology (numlayers) [s]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rheologymodel', 'array describing whether we adopt a Maxwell (0), Burgers (1) or EBM (2) rheology (default: 0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'density', 'array describing each layer\'s density (numlayers) [kg/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'issolid', 'array describing whether the layer is solid or liquid (default: 1) (numlayers)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'earth_density', 'mantle density [kg/m^3]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - materials.litho Class'
        return s

## ------------------------------------------------------
## materials.damageice
## ------------------------------------------------------
@class_registry.register_class
class damageice(class_registry.manage_state):
    '''
    materials.damageice Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.rho_ice = 917.
        self.rho_water = 1023.
        self.rho_freshwater = 1000.
        self.mu_water = 0.001787
        self.heatcapacity = 2093.
        self.latentheat = 3.34 * pow(10, 5)
        self.thermalconductivity = 2.4
        self.temperateiceconductivity = 0.24
        self.effectiveconductivity_averaging = 1
        self.meltingpoint = 273.15
        self.beta = 9.8 * pow(10, -8)
        self.mixed_layer_capacity = 3974.
        self.thermal_exchange_velocity = 1.00e-4
        self.rheology_B = np.nan
        self.rheology_n = np.nan
        self.rheology_law = 'Paterson'
        self.earth_density = 5512.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Materials (damage ice):\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'rho_ice', 'ice density [kg/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rho_water', 'water density [kg/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rho_freshwater', 'fresh water density [kg/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'mu_water', 'water viscosity [N s/m^2]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'heatcapacity', 'heat capacity [J/kg/K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thermalconductivity', 'ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'temperateiceconductivity', 'temperate ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'effectiveconductivity_averaging', 'computation of effectiveconductivity: (0) arithmetic mean, (1) harmonic mean, (2) geometric mean (default)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'meltingpoint', 'melting point of ice at 1atm in K'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'latentheat', 'latent heat of fusion [J/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'beta', 'rate of change of melting point with pressure [K/Pa]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'mixed_layer_capacity', 'mixed layer capacity [W/ kg/K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thermal_exchange_velocity', 'thermal exchange velocity [m/s]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rheology_B', 'flow law parameter [Pa s^(1/n)]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rheology_n', 'Glen\'s flow law exponent'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rheology_law', 'law for the temperature dependance of the rheology: \'None\', \'BuddJacka\', \'Cuffey\', \'CuffeyTemperate\', \'Paterson\', \'Arrhenius\' or \'LliboutryDuval\''))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'earth_density', 'Mantle density [kg m^-3]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - materials.damageice Class'
        return s

## ------------------------------------------------------
## materials.enhancedice
## ------------------------------------------------------
@class_registry.register_class
class enhancedice(class_registry.manage_state):
    '''
    materials.enhancedice Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.rho_ice = 917.
        self.rho_water = 1023.
        self.rho_freshwater = 1000.
        self.mu_water = 0.001787
        self.heatcapacity = 2093.
        self.latentheat = 3.34 * pow(10, 5)
        self.thermalconductivity = 2.4
        self.temperateiceconductivity = 0.24
        self.effectiveconductivity_averaging = 1
        self.meltingpoint = 273.15
        self.beta = 9.8 * pow(10, -8)
        self.mixed_layer_capacity = 3974.
        self.thermal_exchange_velocity = 1.00 * pow(10, -4)
        self.rheology_E = np.nan
        self.rheology_B = np.nan
        self.rheology_n = np.nan
        self.rheology_law = 'Paterson'
        self.earth_density = 5512.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Materials (enhanced ice):\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'rho_ice', 'ice density [kg/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rho_water', 'water density [kg/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rho_freshwater', 'fresh water density [kg/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'mu_water', 'water viscosity [N s/m^2]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'heatcapacity', 'heat capacity [J/kg/K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thermalconductivity', 'ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'temperateiceconductivity', 'temperate ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'effectiveconductivity_averaging', 'computation of effectiveconductivity: (0) arithmetic mean, (1) harmonic mean, (2) geometric mean (default)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'meltingpoint', 'melting point of ice at 1atm in K'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'latentheat', 'latent heat of fusion [J/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'beta', 'rate of change of melting point with pressure [K/Pa]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'mixed_layer_capacity', 'mixed layer capacity [W/kg/K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thermal_exchange_velocity', 'thermal exchange velocity [m/s]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rheology_E', 'enhancement factor'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rheology_B', 'flow law parameter [Pa s^(1/n)]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rheology_n', 'Glen\'s flow law exponent'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rheology_law', 'law for the temperature dependance of the rheology: \'None\', \'BuddJacka\', \'Cuffey\', \'CuffeyTemperate\', \'Paterson\', \'Arrhenius\' or \'LliboutryDuval\''))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'earth_density', 'Mantle density [kg/m^-3]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - materials.enhancedice Class'
        return s

## ------------------------------------------------------
## materials.estar
## ------------------------------------------------------
@class_registry.register_class
class estar(class_registry.manage_state):
    '''
    materials.estar Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.rho_ice = 917.
        self.rho_water = 1023.
        self.rho_freshwater = 1000.
        self.mu_water = 0.001787
        self.heatcapacity = 2093.
        self.latentheat = 3.34 * pow(10, 5)
        self.thermalconductivity = 2.4
        self.temperateiceconductivity = 0.24
        self.effectiveconductivity_averaging = 1
        self.meltingpoint = 273.15
        self.beta = 9.8 * pow(10, -8)
        self.mixed_layer_capacity = 3974.
        self.thermal_exchange_velocity = 1.00 * pow(10, -4)
        self.rheology_B = np.nan
        self.rheology_Ec = np.nan
        self.rheology_Es = np.nan
        self.rheology_law = 'Paterson'
        self.earth_density = 5512.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Materials (estar):\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'rho_ice', 'ice density [kg/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rho_water', 'ocean water density [kg/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rho_freshwater', 'fresh water density [kg/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'mu_water', 'water viscosity [N s/m^2]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'heatcapacity', 'heat capacity [J/kg/K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thermalconductivity', ['ice thermal conductivity [W/m/K]']))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'temperateiceconductivity', 'temperate ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, "effectiveconductivity_averaging", "computation of effectiveconductivity: (0) arithmetic mean, (1) harmonic mean, (2) geometric mean (default)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'meltingpoint', 'melting point of ice at 1atm in K'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'latentheat', 'latent heat of fusion [J/kg]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'beta', 'rate of change of melting point with pressure [K/Pa]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'mixed_layer_capacity', 'mixed layer capacity [W/kg/K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thermal_exchange_velocity', 'thermal exchange velocity [m/s]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rheology_B', 'flow law parameter [Pa s^(1/3)]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rheology_Ec', 'compressive enhancement factor'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rheology_Es', 'shear enhancement factor'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rheology_law', ['law for the temperature dependance of the rheology: \'None\', \'BuddJacka\', \'Cuffey\', \'CuffeyTemperate\', \'Paterson\', \'Arrhenius\' or \'LliboutryDuval\'']))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'earth_density', 'Mantle density [kg/m^-3]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - materials.estar Class'
        return s
