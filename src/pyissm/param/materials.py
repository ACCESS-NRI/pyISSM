import numpy as np
from . import param_utils
from . import class_registry
from .. import execute

## ------------------------------------------------------
## materials.ice
## ------------------------------------------------------
@class_registry.register_class
class ice(class_registry.manage_state):
    """
    Ice materials parameters class for ISSM.

    This class defines the default physical parameters for ice used in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    rho_ice : float, default=917.
        Ice density [kg/m^3]
    rho_water : float, default=1023.
        Ocean water density [kg/m^3]
    rho_freshwater : float, default=1000.
        Fresh water density [kg/m^3]
    mu_water : float, default=0.001787
        Water viscosity [N s/m^2]
    heatcapacity : float, default=2093.
        Heat capacity [J/kg/K]
    latentheat : float, default=3.34e5
        Latent heat of fusion [J/m^3]
    thermalconductivity : float, default=2.4
        Ice thermal conductivity [W/m/K]
    temperateiceconductivity : float, default=0.24
        Temperate ice thermal conductivity [W/m/K]
    effectiveconductivity_averaging : int, default=1
        Computation of effective conductivity: (0) arithmetic mean, (1) harmonic mean, (2) geometric mean
    meltingpoint : float, default=273.15
        Melting point of ice at 1 atm [K]
    beta : float, default=9.8e-8
        Rate of change of melting point with pressure [K/Pa]
    mixed_layer_capacity : float, default=3974.
        Mixed layer capacity [W/kg/K]
    thermal_exchange_velocity : float, default=1.00e-4
        Thermal exchange velocity [m/s]
    rheology_law : str, default='Paterson'
        Law for the temperature dependence of the rheology: 'None', 'BuddJacka', 'Cuffey', 'CuffeyTemperate', 'Paterson', 'Arrhenius', 'LliboutryDuval', 'NyeCO2', or 'NyeH2O'
    rheology_B : float, default=2.1e8
        Flow law parameter [Pa s^(1/n)]
    rheology_n : float, default=3.
        Glen's flow law exponent
    earth_density : float, default=5512.
        Mantle density [kg/m^3]

    Methods
    -------
    __init__(self, other=None)
        Initializes the default ice material parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the ice material parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.materials = pyissm.param.materials.ice()
    """

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

        s += '{}\n'.format(param_utils.fielddisplay(self, 'rho_ice', 'ice density [kg/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rho_water', 'ocean water density [kg/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rho_freshwater', 'fresh water density [kg/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'mu_water', 'water viscosity [N s/m^2]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'heatcapacity', 'heat capacity [J/kg/K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'thermalconductivity', 'ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'temperateiceconductivity', 'temperate ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'meltingpoint', 'melting point of ice at 1atm in K'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'latentheat', 'latent heat of fusion [J/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'beta', 'rate of change of melting point with pressure [K/Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'mixed_layer_capacity', 'mixed layer capacity [W/kg/K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'thermal_exchange_velocity', 'thermal exchange velocity [m/s]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rheology_B', 'flow law parameter [Pa s^(1/n)]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rheology_n', 'Glen\'s flow law exponent'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rheology_law', 'law for the temperature dependance of the rheology: \'None\', \'BuddJacka\', \'Cuffey\', \'CuffeyTemperate\', \'Paterson\', \'Arrhenius\', \'LliboutryDuval\', \'NyeCO2\', or \'NyeH2O\''))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - materials.ice Class'
        return s
    
    # Marshall method for saving the mask parameters
    def marshall_class(self, prefix, md, fid):
        """
        Marshall the mask parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.

        Returns
        -------
        None
        """

        ## Write headers to file
        # NOTE: data types must match the expected types in the ISSM code. These come from naturetointeger in $ISSM_DIR/src/m/materials.py
        execute.WriteData(fid, prefix, data = 3, name = 'md.materials.nature', format = 'IntMat', mattype = 3)
        execute.WriteData(fid, prefix, data = 5, name = 'md.materials.type', format = 'Integer')

        ## Write all 'Double' fields to file
        fieldnames = ['rho_ice', 'rho_water', 'rho_freshwater', 'mu_water', 'heatcapacity',
                      'latentheat', 'thermalconductivity', 'temperateiceconductivity', 'meltingpoint', 'beta',
                      'mixed_layer_capacity', 'thermal_exchange_velocity', 'earth_density']
        for fieldname in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = fieldname, format = 'Double')

        ## Write all other fields to file
        execute.WriteData(fid, prefix, obj = self, fieldname = 'effectiveconductivity_averaging', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rheology_B', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rheology_n', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, data = self.rheology_law, name = 'md.materials.rheology_law', format = 'String')

## ------------------------------------------------------
## materials.hydro
## ------------------------------------------------------
@class_registry.register_class
class hydro(class_registry.manage_state):
    """
    Hydro materials parameters class for ISSM.

    This class defines the default physical parameters for hydro (hydrology) used in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    rho_ice : float, default=917.
        Ice density [kg/m^3]
    rho_water : float, default=1023.
        Ocean water density [kg/m^3]
    rho_freshwater : float, default=1000.
        Fresh water density [kg/m^3]
    earth_density : float, default=5512.
        Mantle density [kg/m^3]

    Methods
    -------
    __init__(self, other=None)
        Initializes the default hydro material parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the hydro material parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.materials = pyissm.param.materials.hydro()
    """

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

        s += '{}\n'.format(param_utils.fielddisplay(self, 'rho_ice', 'ice density [kg/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rho_water', 'ocean water density [kg/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rho_freshwater', 'fresh water density [kg/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'earth_density', 'mantle density [kg/m^3]'))
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
    """
    Lithosphere materials parameters class for ISSM.

    This class defines the default physical parameters for the lithosphere used in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    numlayers : int, default=2
        Number of layers in the lithosphere model.
    radius : list of float, default=[1e3, 6278e3, 6378e3]
        Radii for each interface (numlayers + 1) [m].
    viscosity : list of float, default=[1e21, 1e40]
        Viscosity for each layer (numlayers) [Pa.s].
    lame_mu : list of float, default=[1.45e11, 6.7e10]
        Shear modulus for each layer (numlayers) [Pa].
    lame_lambda : list of float, default=[1.45e11, 6.7e10]
        Lame lambda parameter for each layer (numlayers) [Pa].
    burgers_viscosity : list of float, default=[np.nan, np.nan]
        Transient viscosity for Burgers rheologies (numlayers) [Pa.s].
    burgers_mu : list of float, default=[np.nan, np.nan]
        Transient shear modulus for Burgers rheologies (numlayers) [Pa].
    ebm_alpha : list of float, default=[np.nan, np.nan]
        Exponent parameter for EBM rheology (numlayers).
    ebm_delta : list of float, default=[np.nan, np.nan]
        Amplitude of transient relaxation for EBM rheology (numlayers).
    ebm_taul : list of float, default=[np.nan, np.nan]
        Starting period for transient relaxation for EBM rheology (numlayers) [s].
    ebm_tauh : list of float, default=[np.nan, np.nan]
        End period for transient relaxation for Burgers rheology (numlayers) [s].
    rheologymodel : list of int, default=[0, 0]
        Rheology model for each layer: Maxwell (0), Burgers (1), or EBM (2).
    density : list of float, default=[5.51e3, 5.50e3]
        Density for each layer (numlayers) [kg/m^3].
    issolid : list of int, default=[1, 1]
        Whether each layer is solid (1) or liquid (0) (numlayers).
    earth_density : float, default=5512.
        Mantle density [kg/m^3].

    Methods
    -------
    __init__(self, other=None)
        Initializes the default lithosphere material parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the lithosphere material parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.materials = pyissm.param.materials.litho()
    """

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

        s += '{}\n'.format(param_utils.fielddisplay(self, 'numlayers', 'number of layers (default: 2)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'radius', 'array describing the radius for each interface (numlayers + 1) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'viscosity', 'array describing each layer\'s viscosity (numlayers) [Pa.s]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'lame_lambda', 'array describing the lame lambda parameter (numlayers) [Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'lame_mu', 'array describing the shear modulus for each layers (numlayers) [Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'burgers_viscosity', 'array describing each layer\'s transient viscosity, only for Burgers rheologies  (numlayers) [Pa.s]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'burgers_mu', 'array describing each layer\'s transient shear modulus, only for Burgers rheologies  (numlayers) [Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ebm_alpha', 'array describing each layer\'s exponent parameter controlling the shape of shear modulus curve between taul and tauh, only for EBM rheology (numlayers)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ebm_delta', 'array describing each layer\'s amplitude of the transient relaxation (ratio between elastic rigity to pre-maxwell relaxation rigity), only for EBM rheology (numlayers)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ebm_taul', 'array describing each layer\'s starting period for transient relaxation, only for EBM rheology  (numlayers) [s]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ebm_tauh', 'array describing each layer''s array describing each layer\'s end period for transient relaxation, only for Burgers rheology (numlayers) [s]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rheologymodel', 'array describing whether we adopt a Maxwell (0), Burgers (1) or EBM (2) rheology (default: 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'density', 'array describing each layer\'s density (numlayers) [kg/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'issolid', 'array describing whether the layer is solid or liquid (default: 1) (numlayers)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'earth_density', 'mantle density [kg/m^3]'))
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
    """
    Damage ice materials parameters class for ISSM.

    This class defines the default physical parameters for damage ice used in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    rho_ice : float, default=917.
        Ice density [kg/m^3]
    rho_water : float, default=1023.
        Ocean water density [kg/m^3]
    rho_freshwater : float, default=1000.
        Fresh water density [kg/m^3]
    mu_water : float, default=0.001787
        Water viscosity [N s/m^2]
    heatcapacity : float, default=2093.
        Heat capacity [J/kg/K]
    latentheat : float, default=3.34e5
        Latent heat of fusion [J/m^3]
    thermalconductivity : float, default=2.4
        Ice thermal conductivity [W/m/K]
    temperateiceconductivity : float, default=0.24
        Temperate ice thermal conductivity [W/m/K]
    effectiveconductivity_averaging : int, default=1
        Computation of effective conductivity: (0) arithmetic mean, (1) harmonic mean, (2) geometric mean
    meltingpoint : float, default=273.15
        Melting point of ice at 1 atm [K]
    beta : float, default=9.8e-8
        Rate of change of melting point with pressure [K/Pa]
    mixed_layer_capacity : float, default=3974.
        Mixed layer capacity [W/kg/K]
    thermal_exchange_velocity : float, default=1.00e-4
        Thermal exchange velocity [m/s]
    rheology_B : float, default=np.nan
        Flow law parameter [Pa s^(1/n)]
    rheology_n : float, default=np.nan
        Glen's flow law exponent
    rheology_law : str, default='Paterson'
        Law for the temperature dependence of the rheology: 'None', 'BuddJacka', 'Cuffey', 'CuffeyTemperate', 'Paterson', 'Arrhenius', 'LliboutryDuval', 'NyeCO2', or 'NyeH2O'
    earth_density : float, default=5512.
        Mantle density [kg/m^3]

    Methods
    -------
    __init__(self, other=None)
        Initializes the default damage ice material parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the damage ice material parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.materials = pyissm.param.materials.damageice()
    """

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

        s += '{}\n'.format(param_utils.fielddisplay(self, 'rho_ice', 'ice density [kg/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rho_water', 'water density [kg/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rho_freshwater', 'fresh water density [kg/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'mu_water', 'water viscosity [N s/m^2]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'heatcapacity', 'heat capacity [J/kg/K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'thermalconductivity', 'ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'temperateiceconductivity', 'temperate ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'effectiveconductivity_averaging', 'computation of effectiveconductivity: (0) arithmetic mean, (1) harmonic mean, (2) geometric mean (default)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'meltingpoint', 'melting point of ice at 1atm in K'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'latentheat', 'latent heat of fusion [J/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'beta', 'rate of change of melting point with pressure [K/Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'mixed_layer_capacity', 'mixed layer capacity [W/ kg/K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'thermal_exchange_velocity', 'thermal exchange velocity [m/s]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rheology_B', 'flow law parameter [Pa s^(1/n)]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rheology_n', 'Glen\'s flow law exponent'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rheology_law', 'law for the temperature dependance of the rheology: \'None\', \'BuddJacka\', \'Cuffey\', \'CuffeyTemperate\', \'Paterson\', \'Arrhenius\' or \'LliboutryDuval\''))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'earth_density', 'Mantle density [kg m^-3]'))
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
    """
    Enhanced ice materials parameters class for ISSM.

    This class defines the default physical parameters for enhanced ice used in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    rho_ice : float, default=917.
        Ice density [kg/m^3]
    rho_water : float, default=1023.
        Ocean water density [kg/m^3]
    rho_freshwater : float, default=1000.
        Fresh water density [kg/m^3]
    mu_water : float, default=0.001787
        Water viscosity [N s/m^2]
    heatcapacity : float, default=2093.
        Heat capacity [J/kg/K]
    latentheat : float, default=3.34e5
        Latent heat of fusion [J/m^3]
    thermalconductivity : float, default=2.4
        Ice thermal conductivity [W/m/K]
    temperateiceconductivity : float, default=0.24
        Temperate ice thermal conductivity [W/m/K]
    effectiveconductivity_averaging : int, default=1
        Computation of effective conductivity: (0) arithmetic mean, (1) harmonic mean, (2) geometric mean
    meltingpoint : float, default=273.15
        Melting point of ice at 1 atm [K]
    beta : float, default=9.8e-8
        Rate of change of melting point with pressure [K/Pa]
    mixed_layer_capacity : float, default=3974.
        Mixed layer capacity [W/kg/K]
    thermal_exchange_velocity : float, default=1.00e-4
        Thermal exchange velocity [m/s]
    rheology_E : float, default=np.nan
        Enhancement factor
    rheology_B : float, default=np.nan
        Flow law parameter [Pa s^(1/n)]
    rheology_n : float, default=np.nan
        Glen's flow law exponent
    rheology_law : str, default='Paterson'
        Law for the temperature dependence of the rheology: 'None', 'BuddJacka', 'Cuffey', 'CuffeyTemperate', 'Paterson', 'Arrhenius', or 'LliboutryDuval'
    earth_density : float, default=5512.
        Mantle density [kg/m^3]

    Methods
    -------
    __init__(self, other=None)
        Initializes the default enhanced ice material parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the enhanced ice material parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.materials = pyissm.param.materials.enhancedice()
    """

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

        s += '{}\n'.format(param_utils.fielddisplay(self, 'rho_ice', 'ice density [kg/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rho_water', 'water density [kg/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rho_freshwater', 'fresh water density [kg/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'mu_water', 'water viscosity [N s/m^2]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'heatcapacity', 'heat capacity [J/kg/K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'thermalconductivity', 'ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'temperateiceconductivity', 'temperate ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'effectiveconductivity_averaging', 'computation of effectiveconductivity: (0) arithmetic mean, (1) harmonic mean, (2) geometric mean (default)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'meltingpoint', 'melting point of ice at 1atm in K'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'latentheat', 'latent heat of fusion [J/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'beta', 'rate of change of melting point with pressure [K/Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'mixed_layer_capacity', 'mixed layer capacity [W/kg/K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'thermal_exchange_velocity', 'thermal exchange velocity [m/s]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rheology_E', 'enhancement factor'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rheology_B', 'flow law parameter [Pa s^(1/n)]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rheology_n', 'Glen\'s flow law exponent'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rheology_law', 'law for the temperature dependance of the rheology: \'None\', \'BuddJacka\', \'Cuffey\', \'CuffeyTemperate\', \'Paterson\', \'Arrhenius\' or \'LliboutryDuval\''))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'earth_density', 'Mantle density [kg/m^-3]'))
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
    """
    E* (estar) ice materials parameters class for ISSM.

    This class defines the default physical parameters for E* (estar) ice used in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    rho_ice : float, default=917.
        Ice density [kg/m^3]
    rho_water : float, default=1023.
        Ocean water density [kg/m^3]
    rho_freshwater : float, default=1000.
        Fresh water density [kg/m^3]
    mu_water : float, default=0.001787
        Water viscosity [N s/m^2]
    heatcapacity : float, default=2093.
        Heat capacity [J/kg/K]
    latentheat : float, default=3.34e5
        Latent heat of fusion [J/m^3]
    thermalconductivity : float, default=2.4
        Ice thermal conductivity [W/m/K]
    temperateiceconductivity : float, default=0.24
        Temperate ice thermal conductivity [W/m/K]
    effectiveconductivity_averaging : int, default=1
        Computation of effective conductivity: (0) arithmetic mean, (1) harmonic mean, (2) geometric mean
    meltingpoint : float, default=273.15
        Melting point of ice at 1 atm [K]
    beta : float, default=9.8e-8
        Rate of change of melting point with pressure [K/Pa]
    mixed_layer_capacity : float, default=3974.
        Mixed layer capacity [W/kg/K]
    thermal_exchange_velocity : float, default=1.00e-4
        Thermal exchange velocity [m/s]
    rheology_B : ndarray, default=np.nan
        Flow law parameter [Pa s^(1/3)]
    rheology_Ec : ndarray, default=np.nan
        Compressive enhancement factor
    rheology_Es : ndarray, default=np.nan
        Shear enhancement factor
    rheology_law : str, default='Paterson'
        Law for the temperature dependence of the rheology: 'None', 'BuddJacka', 'Cuffey', 'CuffeyTemperate', 'Paterson', 'Arrhenius', or 'LliboutryDuval'
    earth_density : float, default=5512.
        Mantle density [kg/m^3]

    Methods
    -------
    __init__(self, other=None)
        Initializes the default E* ice material parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the E* ice material parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.materials = pyissm.param.materials.estar()
    """

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

        s += '{}\n'.format(param_utils.fielddisplay(self, 'rho_ice', 'ice density [kg/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rho_water', 'ocean water density [kg/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rho_freshwater', 'fresh water density [kg/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'mu_water', 'water viscosity [N s/m^2]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'heatcapacity', 'heat capacity [J/kg/K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'thermalconductivity', ['ice thermal conductivity [W/m/K]']))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'temperateiceconductivity', 'temperate ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, "effectiveconductivity_averaging", "computation of effectiveconductivity: (0) arithmetic mean, (1) harmonic mean, (2) geometric mean (default)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'meltingpoint', 'melting point of ice at 1atm in K'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'latentheat', 'latent heat of fusion [J/kg]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'beta', 'rate of change of melting point with pressure [K/Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'mixed_layer_capacity', 'mixed layer capacity [W/kg/K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'thermal_exchange_velocity', 'thermal exchange velocity [m/s]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rheology_B', 'flow law parameter [Pa s^(1/3)]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rheology_Ec', 'compressive enhancement factor'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rheology_Es', 'shear enhancement factor'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rheology_law', ['law for the temperature dependance of the rheology: \'None\', \'BuddJacka\', \'Cuffey\', \'CuffeyTemperate\', \'Paterson\', \'Arrhenius\' or \'LliboutryDuval\'']))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'earth_density', 'Mantle density [kg/m^-3]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - materials.estar Class'
        return s
