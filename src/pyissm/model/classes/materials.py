import numpy as np
import warnings
from dataclasses import dataclass
from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm.model import execute, mesh

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
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.materials = pyissm.model.classes.materials.ice()
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

        s += '{}\n'.format(class_utils.fielddisplay(self, 'rho_ice', 'ice density [kg/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rho_water', 'ocean water density [kg/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rho_freshwater', 'fresh water density [kg/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'mu_water', 'water viscosity [N s/m^2]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'heatcapacity', 'heat capacity [J/kg/K]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'thermalconductivity', 'ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'temperateiceconductivity', 'temperate ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'meltingpoint', 'melting point of ice at 1atm in K'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'latentheat', 'latent heat of fusion [J/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'beta', 'rate of change of melting point with pressure [K/Pa]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'mixed_layer_capacity', 'mixed layer capacity [W/kg/K]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'thermal_exchange_velocity', 'thermal exchange velocity [m/s]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rheology_B', 'flow law parameter [Pa s^(1/n)]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rheology_n', 'Glen\'s flow law exponent'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rheology_law', 'law for the temperature dependance of the rheology: \'None\', \'BuddJacka\', \'Cuffey\', \'CuffeyTemperate\', \'Paterson\', \'Arrhenius\', \'LliboutryDuval\', \'NyeCO2\', or \'NyeH2O\''))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - materials.ice Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [materials.ice] fields to 3D
        """
        self.rheology_B = mesh.project_3d(md, vector = self.rheology_B, type = 'node')
        self.rheology_n = mesh.project_3d(md, vector = self.rheology_n, type = 'element')

        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        if solution == 'TransientSolution' and md.transient.isslc:
            class_utils.check_field(md, fieldname = 'materials.earth_density', scalar = True, gt = 0)
        else:
            class_utils.check_field(md, fieldname = 'materials.rho_ice', gt = 0)
            class_utils.check_field(md, fieldname = 'materials.rho_water', gt = 0)
            class_utils.check_field(md, fieldname = 'materials.rho_freshwater', gt = 0)
            class_utils.check_field(md, fieldname = 'materials.mu_water', gt = 0)
            class_utils.check_field(md, fieldname = 'materials.rheology_B', gt = 0, size = 'universal', allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = 'materials.rheology_n', gt = 0, size = 'universal', allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = 'materials.rheology_law', values = ['None', 'BuddJacka', 'Cuffey', 'CuffeyTemperate', 'Paterson', 'Arrhenius', 'LliboutryDuval', 'NyeCO2', 'NyeH2O'])
            class_utils.check_field(md, fieldname = 'materials.effectiveconductivity_averaging', scalar = True, values = [0, 1, 2])

        return md
    
    # Marshall method for saving the materials.ice parameters
    def marshall_class(self, fid, prefix, md = None, write_type = True):
        """
        Marshall [materials.ice] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
            ISSM model object needed in some cases.

        Returns
        -------
        None
        """

        ## Write headers to file
        if write_type:
            execute.WriteData(fid, prefix, data = _material_registry['ice'].type, name = 'md.materials.type', format = 'Integer')

        ## Write Double fields
        fieldnames = ['rho_ice', 'rho_water', 'rho_freshwater', 'mu_water', 'heatcapacity',
                      'latentheat', 'thermalconductivity', 'temperateiceconductivity', 'meltingpoint', 'beta',
                      'mixed_layer_capacity', 'thermal_exchange_velocity', 'earth_density']
        for fieldname in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = fieldname, format = 'Double')

        ## Write conditional fields
        if (
            np.size(self.rheology_B) == 1
            or (
                self.rheology_B.ndim == 1
                and self.rheology_B.shape[0] in (
                    md.mesh.numberofvertices,
                    md.mesh.numberofvertices + 1
                )
            )
            or (
                self.rheology_B.ndim == 2
                and (self.rheology_B.shape[0] == md.mesh.numberofvertices or 
                self.rheology_B.shape[0] == md.mesh.numberofvertices + 1 or 
                (self.rheology_B.shape[0] == md.mesh.numberofelements and self.rheology_B.shape[1] > 1)
                )
            )
        ):
            mattype = 1
            tsl = md.mesh.numberofvertices
        else:
            mattype = 2
            tsl = md.mesh.numberofelements
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rheology_B', format = 'DoubleMat', mattype = mattype, timeserieslength = tsl + 1, yts = md.constants.yts)

        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'effectiveconductivity_averaging', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rheology_n', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rheology_law', format = 'String')

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
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.materials = pyissm.model.classes.materials.hydro()
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

        s += '{}\n'.format(class_utils.fielddisplay(self, 'rho_ice', 'ice density [kg/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rho_water', 'ocean water density [kg/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rho_freshwater', 'fresh water density [kg/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'earth_density', 'mantle density [kg/m^3]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - materials.hydro Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [materials.hydro] fields to 3D
        """
        warnings.warn('pyissm.model.classes.materials.hydro._extrude: 3D extrusion not implemented for materials.hydro. Returning unchanged (2D) materials fields.')

        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        class_utils.check_field(md, fieldname = 'materials.rho_ice', gt = 0)
        class_utils.check_field(md, fieldname = 'materials.rho_water', gt = 0)
        class_utils.check_field(md, fieldname = 'materials.earth_density', scalar = True, gt = 0)
        class_utils.check_field(md, fieldname = 'materials.rho_freshwater', gt = 0)

        return md
    
    # Marshall method for saving the materials.hydro parameters
    def marshall_class(self, fid, prefix, md = None, write_type = True):
        """
        Marshall [materials.hydro] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
            ISSM model object needed in some cases.

        Returns
        -------
        None
        """

        ## Write headers to file
        if write_type:
            execute.WriteData(fid, prefix, data = _material_registry['hydro'].type, name = 'md.materials.type', format = 'Integer')

        ## Write fields (all consistent format)
        fieldnames = ['rho_ice', 'rho_water', 'rho_freshwater', 'earth_density']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Double')

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
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.materials = pyissm.model.classes.materials.litho()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.numlayers = 2
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

        s += '{}\n'.format(class_utils.fielddisplay(self, 'numlayers', 'number of layers (default: 2)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'radius', 'array describing the radius for each interface (numlayers + 1) [m]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'viscosity', 'array describing each layer\'s viscosity (numlayers) [Pa.s]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'lame_lambda', 'array describing the lame lambda parameter (numlayers) [Pa]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'lame_mu', 'array describing the shear modulus for each layers (numlayers) [Pa]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'burgers_viscosity', 'array describing each layer\'s transient viscosity, only for Burgers rheologies  (numlayers) [Pa.s]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'burgers_mu', 'array describing each layer\'s transient shear modulus, only for Burgers rheologies  (numlayers) [Pa]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'ebm_alpha', 'array describing each layer\'s exponent parameter controlling the shape of shear modulus curve between taul and tauh, only for EBM rheology (numlayers)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'ebm_delta', 'array describing each layer\'s amplitude of the transient relaxation (ratio between elastic rigity to pre-maxwell relaxation rigity), only for EBM rheology (numlayers)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'ebm_taul', 'array describing each layer\'s starting period for transient relaxation, only for EBM rheology  (numlayers) [s]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'ebm_tauh', 'array describing each layer''s array describing each layer\'s end period for transient relaxation, only for Burgers rheology (numlayers) [s]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rheologymodel', 'array describing whether we adopt a Maxwell (0), Burgers (1) or EBM (2) rheology (default: 0)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'density', 'array describing each layer\'s density (numlayers) [kg/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'issolid', 'array describing whether the layer is solid or liquid (default: 1) (numlayers)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'earth_density', 'mantle density [kg/m^3]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - materials.litho Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [materials.litho] fields to 3D
        """
        warnings.warn('pyissm.model.classes.materials.litho._extrude: 3D extrusion not implemented for materials.litho. Returning unchanged (2D) materials fields.')

        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        # Early return if not LoveAnalysis
        if 'LoveAnalysis' not in analyses:
            return md
        
        class_utils.check_field(md, fieldname = 'materials.numlayers', scalar = True, gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'materials.radius', size = (md.materials.numlayers + 1, 1), gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'materials.lame_mu', size = (md.materials.numlayers, 1), ge = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'materials.lame_lambda', size = (md.materials.numlayers, 1), ge = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'materials.issolid', size = (md.materials.numlayers, 1), ge = 0, lt = 2, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'materials.density', size = (md.materials.numlayers, 1), gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'materials.viscosity', size = (md.materials.numlayers, 1), ge = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'materials.rheologymodel', size = (md.materials.numlayers, 1), ge = 0, le = 2, allow_nan = False, allow_inf = False)

        if np.any(self.rheologymodel == 1):
            class_utils.check_field(md, fieldname = 'materials.burgers_viscosity', size = (md.materials.numlayers, 1), ge = 0, allow_inf = False)
            class_utils.check_field(md, fieldname = 'materials.burgers_mu', size = (md.materials.numlayers, 1), ge = 0, allow_inf = False)

        if np.any(self.rheologymodel == 2):
            class_utils.check_field(md, fieldname = 'materials.ebm_alpha', size = (md.materials.numlayers, 1), ge = 0, allow_inf = False)
            class_utils.check_field(md, fieldname = 'materials.ebm_delta', size = (md.materials.numlayers, 1), ge = 0, allow_inf = False)
            class_utils.check_field(md, fieldname = 'materials.ebm_taul', size = (md.materials.numlayers, 1), ge = 0, allow_inf = False)
            class_utils.check_field(md, fieldname = 'materials.ebm_tauh', size = (md.materials.numlayers, 1), ge = 0, allow_inf = False)

        for i in range(md.materials.numlayers):
            if md.materials.rheologymodel[i] == 1 and (np.isnan(md.materials.burgers_viscosity[i] or np.isnan(md.materials.burgers_mu[i]))):
                raise RuntimeError('pyissm.model.classes.materials.litho.check_consistency: Litho burgers_viscosity or burgers_mu has NaN values, inconsistent with rheologymodel choice')

            if md.materials.rheologymodel[i] == 2 and (np.isnan(md.materials.ebm_alpha[i]) or np.isnan(md.materials.ebm_delta[i]) or np.isnan(md.materials.ebm_taul[i]) or np.isnan(md.materials.ebm_tauh[i])):
                raise RuntimeError('pyissm.model.classes.materials.litho.check_consistency: Litho ebm_alpha, ebm_delta, ebm_taul or ebm_tauh has NaN values, inconsistent with rheologymodel choice')
        if md.materials.issolid[0] == 0 or md.materials.lame_mu[0] == 0:
            raise RuntimeError('First layer must be solid (issolid[0] > 0 AND lame_mu[0] > 0). Add a weak inner core if necessary.')
        ind = np.where(md.materials.issolid == 0)[0]

        return md
    
    # Marshall method for saving the materials.litho parameters
    def marshall_class(self, fid, prefix, md = None, write_type = True):
        """
        Marshall [materials.litho] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
            ISSM model object needed in some cases.

        Returns
        -------
        None
        """

        ## Compute earth density compatible with our layer density distribution
        earth_density = 0
        for i in range(self.numlayers):
            earth_density = earth_density + (pow(self.radius[i + 1], 3) - pow(self.radius[i], 3)) * self.density[i]
        earth_density = earth_density / pow(self.radius[self.numlayers], 3)
        self.earth_density = earth_density

        ## Write headers to file
        if write_type:
            execute.WriteData(fid, prefix, data = _material_registry['litho'].type, name = 'md.materials.type', format = 'Integer')

        ## Write DoubleMat fields
        fieldnames = ['radius', 'lame_mu', 'lame_lambda', 'issolid', 'density', 'viscosity', 'rheologymodel',
                      'burgers_viscosity', 'burgers_mu', 'ebm_alpha', 'ebm_delta', 'ebm_taul', 'ebm_tauh']
        for field in fieldnames:
                execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'DoubleMat', mattype = 3)
        
        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'numlayers', format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.materials.earth_density', data = self.earth_density, format = 'Double')

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
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.materials = pyissm.model.classes.materials.damageice()
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

        s += '{}\n'.format(class_utils.fielddisplay(self, 'rho_ice', 'ice density [kg/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rho_water', 'water density [kg/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rho_freshwater', 'fresh water density [kg/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'mu_water', 'water viscosity [N s/m^2]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'heatcapacity', 'heat capacity [J/kg/K]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'thermalconductivity', 'ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'temperateiceconductivity', 'temperate ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'effectiveconductivity_averaging', 'computation of effectiveconductivity: (0) arithmetic mean, (1) harmonic mean, (2) geometric mean (default)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'meltingpoint', 'melting point of ice at 1atm in K'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'latentheat', 'latent heat of fusion [J/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'beta', 'rate of change of melting point with pressure [K/Pa]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'mixed_layer_capacity', 'mixed layer capacity [W/ kg/K]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'thermal_exchange_velocity', 'thermal exchange velocity [m/s]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rheology_B', 'flow law parameter [Pa s^(1/n)]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rheology_n', 'Glen\'s flow law exponent'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rheology_law', 'law for the temperature dependance of the rheology: \'None\', \'BuddJacka\', \'Cuffey\', \'CuffeyTemperate\', \'Paterson\', \'Arrhenius\' or \'LliboutryDuval\''))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'earth_density', 'Mantle density [kg m^-3]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - materials.damageice Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [materials.damageice] fields to 3D
        """
        self.rheology_B = mesh.project_3d(md, vector = self.rheology_B, type = 'node')
        self.rheology_n = mesh.project_3d(md, vector = self.rheology_n, type = 'element')

        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        class_utils.check_field(md, fieldname = 'materials.rho_ice', gt = 0)
        class_utils.check_field(md, fieldname = 'materials.rho_water', gt = 0)
        class_utils.check_field(md, fieldname = 'materials.rho_freshwater', gt = 0)
        class_utils.check_field(md, fieldname = 'materials.mu_water', gt = 0)
        class_utils.check_field(md, fieldname = 'materials.rheology_B', size = (md.mesh.numberofvertices, ), gt = 0)
        class_utils.check_field(md, fieldname = 'materials.rheology_n', size = (md.mesh.numberofelements, ), gt = 0)
        class_utils.check_field(md, fieldname = 'materials.rheology_law', values = ['None', 'BuddJacka', 'Cuffey', 'CuffeyTemperate', 'Paterson', 'Arrhenius', 'LliboutryDuval'])
        class_utils.check_field(md, fieldname = 'materials.effectiveconductivity_averaging', scalar = True, values = [0, 1, 2])

        if 'SealevelchangeAnalysis' in analyses:
            class_utils.check_field(md, fieldname = 'materials.earth_density', scalar = True, gt = 0)

        return md
    
    # Marshall method for saving the materials.damageice parameters
    def marshall_class(self, fid, prefix, md = None, write_type = True):
        """
        Marshall [materials.damageice] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
            ISSM model object needed in some cases.

        Returns
        -------
        None
        """

        ## Write headers to file
        if write_type:
            execute.WriteData(fid, prefix, data = _material_registry['damageice'].type, name = 'md.materials.type', format = 'Integer')

        ## Write Double fields
        fieldnames = ['rho_ice', 'rho_water', 'rho_freshwater', 'mu_water', 'heatcapacity', 'thermalconductivity',
                      'temperateiceconductivity', 'meltingpoint', 'latentheat',
                      'beta', 'mixed_layer_capacity', 'thermal_exchange_velocity', 'earth_density']
        for field in fieldnames:
                execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Double')
        
        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'effectiveconductivity_averaging', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rheology_B', format = 'DoubleMat', mattype = 1, timeserieslength =  md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rheology_n', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rheology_law', format = 'String')

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
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.materials = pyissm.model.classes.materials.enhancedice()
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

        s += '{}\n'.format(class_utils.fielddisplay(self, 'rho_ice', 'ice density [kg/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rho_water', 'water density [kg/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rho_freshwater', 'fresh water density [kg/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'mu_water', 'water viscosity [N s/m^2]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'heatcapacity', 'heat capacity [J/kg/K]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'thermalconductivity', 'ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'temperateiceconductivity', 'temperate ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'effectiveconductivity_averaging', 'computation of effectiveconductivity: (0) arithmetic mean, (1) harmonic mean, (2) geometric mean (default)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'meltingpoint', 'melting point of ice at 1atm in K'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'latentheat', 'latent heat of fusion [J/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'beta', 'rate of change of melting point with pressure [K/Pa]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'mixed_layer_capacity', 'mixed layer capacity [W/kg/K]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'thermal_exchange_velocity', 'thermal exchange velocity [m/s]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rheology_E', 'enhancement factor'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rheology_B', 'flow law parameter [Pa s^(1/n)]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rheology_n', 'Glen\'s flow law exponent'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rheology_law', 'law for the temperature dependance of the rheology: \'None\', \'BuddJacka\', \'Cuffey\', \'CuffeyTemperate\', \'Paterson\', \'Arrhenius\' or \'LliboutryDuval\''))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'earth_density', 'Mantle density [kg/m^-3]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - materials.enhancedice Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [materials.enhancedice] fields to 3D
        """
        self.rheology_E = mesh.project_3d(md, vector = self.rheology_E, type = 'node')
        self.rheology_B = mesh.project_3d(md, vector = self.rheology_B, type = 'node')
        self.rheology_n = mesh.project_3d(md, vector = self.rheology_n, type = 'element')

        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        class_utils.check_field(md, fieldname = 'materials.rho_ice', gt = 0)
        class_utils.check_field(md, fieldname = 'materials.rho_water', gt = 0)
        class_utils.check_field(md, fieldname = 'materials.rho_freshwater', gt = 0)
        class_utils.check_field(md, fieldname = 'materials.mu_water', gt = 0)
        class_utils.check_field(md, fieldname = 'materials.rheology_E', timeseries = True, gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'materials.rheology_B', timeseries = True, gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'materials.rheology_n', size = (md.mesh.numberofelements, ), gt = 0)
        class_utils.check_field(md, fieldname = 'materials.rheology_law', values = ['None', 'BuddJacka', 'Cuffey', 'CuffeyTemperate', 'Paterson', 'Arrhenius', 'LliboutryDuval'])
        class_utils.check_field(md, fieldname = 'materials.effectiveconductivity_averaging', scalar = True, values = [0, 1, 2])

        if 'SealevelchangeAnalysis' in analyses:
            class_utils.check_field(md, fieldname = 'materials.earth_density', scalar = True, gt = 0)

        return md
    
    # Marshall method for saving the materials.enhancedice parameters
    def marshall_class(self, fid, prefix, md = None, write_type = True):
        """
        Marshall [materials.enhancedice] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
            ISSM model object needed in some cases.

        Returns
        -------
        None
        """

        ## Write headers to file
        if write_type:
            execute.WriteData(fid, prefix, data = _material_registry['enhancedice'].type, name = 'md.materials.type', format = 'Integer')

        ## Write Double fields
        fieldnames = ['rho_ice', 'rho_water', 'rho_freshwater', 'mu_water', 'heatcapacity', 'thermalconductivity',
                      'temperateiceconductivity', 'meltingpoint', 'latentheat',
                      'beta', 'mixed_layer_capacity', 'thermal_exchange_velocity', 'earth_density']
        for field in fieldnames:
                execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Double')
        
        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'effectiveconductivity_averaging', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rheology_E', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rheology_B', format = 'DoubleMat', mattype = 1, timeserieslength =  md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rheology_n', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rheology_law', format = 'String')

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
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.materials = pyissm.model.classes.materials.estar()
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

        s += '{}\n'.format(class_utils.fielddisplay(self, 'rho_ice', 'ice density [kg/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rho_water', 'ocean water density [kg/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rho_freshwater', 'fresh water density [kg/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'mu_water', 'water viscosity [N s/m^2]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'heatcapacity', 'heat capacity [J/kg/K]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'thermalconductivity', ['ice thermal conductivity [W/m/K]']))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'temperateiceconductivity', 'temperate ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, "effectiveconductivity_averaging", "computation of effectiveconductivity: (0) arithmetic mean, (1) harmonic mean, (2) geometric mean (default)"))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'meltingpoint', 'melting point of ice at 1atm in K'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'latentheat', 'latent heat of fusion [J/kg]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'beta', 'rate of change of melting point with pressure [K/Pa]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'mixed_layer_capacity', 'mixed layer capacity [W/kg/K]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'thermal_exchange_velocity', 'thermal exchange velocity [m/s]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rheology_B', 'flow law parameter [Pa s^(1/3)]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rheology_Ec', 'compressive enhancement factor'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rheology_Es', 'shear enhancement factor'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'rheology_law', ['law for the temperature dependance of the rheology: \'None\', \'BuddJacka\', \'Cuffey\', \'CuffeyTemperate\', \'Paterson\', \'Arrhenius\' or \'LliboutryDuval\'']))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'earth_density', 'Mantle density [kg/m^-3]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - materials.estar Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [materials.estar] fields to 3D
        """
        self.rheology_B = mesh.project_3d(md, vector = self.rheology_B, type = 'node')
        self.rheology_Ec = mesh.project_3d(md, vector = self.rheology_Ec, type = 'node')
        self.rheology_Es = mesh.project_3d(md, vector = self.rheology_Es, type = 'node')

        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        class_utils.check_field(md, fieldname = 'materials.rho_ice', gt = 0)
        class_utils.check_field(md, fieldname = 'materials.rho_water', gt = 0)
        class_utils.check_field(md, fieldname = 'materials.rho_freshwater', gt = 0)
        class_utils.check_field(md, fieldname = 'materials.mu_water', gt = 0)
        class_utils.check_field(md, fieldname = 'materials.rheology_B', size = (md.mesh.numberofvertices, ), gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'materials.rheology_Ec', size = (md.mesh.numberofvertices, ), gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'materials.rheology_Es', size = (md.mesh.numberofvertices, ), gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'materials.rheology_law', values = ['None', 'BuddJacka', 'Cuffey', 'CuffeyTemperate', 'Paterson', 'Arrhenius', 'LliboutryDuval'])
        class_utils.check_field(md, fieldname = 'materials.effectiveconductivity_averaging', scalar = True, values = [0, 1, 2])

        if 'SealevelchangeAnalysis' in analyses:
            class_utils.check_field(md, fieldname = 'materials.earth_density', scalar = True, gt = 0)

        return md

    # Marshall method for saving the materials.estar parameters
    def marshall_class(self, fid, prefix, md = None, write_type = True):
        """
        Marshall [materials.estar] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
            ISSM model object needed in some cases.

        Returns
        -------
        None
        """

        ## Write headers to file
        if write_type:
            execute.WriteData(fid, prefix, data = _material_registry['estar'].type, name = 'md.materials.type', format = 'Integer')

        ## Write Double fields
        fieldnames = ['rho_ice', 'rho_water', 'rho_freshwater', 'mu_water', 'heatcapacity', 'thermalconductivity',
                      'temperateiceconductivity', 'meltingpoint', 'latentheat',
                      'beta', 'mixed_layer_capacity', 'thermal_exchange_velocity', 'earth_density']
        for field in fieldnames:
                execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Double')
        
        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'effectiveconductivity_averaging', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rheology_B', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rheology_Ec', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rheology_Es', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rheology_law', format = 'String')

## ------------------------------------------------------
## materials.composite
## ------------------------------------------------------
@class_registry.register_class
class composite(class_registry.manage_state):

    def __init__(self, *args, **kwargs):
        raise NotImplementedError('pyissm.model.classes.materials.composite(): Not yet implemented.')

## ------------------------------------------------------
## materials._material_specs (internal class only)
## ------------------------------------------------------
@dataclass(frozen = True)
class _material_specs:
    cls: type
    type: int
    nature: int

## ------------------------------------------------------
## materials._material_registry (internal dict only)
## ------------------------------------------------------
## Define registry of all material classes.
## NOTE: data types must match the expected types in the ISSM code. See $ISS_DIR/src/c/shared/io/Marshalling/IoCodeConversions.cpp 
## composite() replaces materials() and uses np.nan for the nature as it's comprised of various other nature values.
_material_registry = {
    'damageice': _material_specs(cls = damageice,
                                type = 1,
                                nature = 1),
    'estar': _material_specs(cls = estar,
                            type = 2,
                            nature = 2),
    'ice': _material_specs(cls = ice,
                          type = 3,
                          nature = 3),
    'enhancedice': _material_specs(cls = enhancedice,
                                  type = 4,
                                  nature = 4),
    'composite': _material_specs(cls = composite,
                               type = 5,
                               nature = np.nan),
    'litho': _material_specs(cls = litho,
                            type = 6,
                            nature = 6),
    'hydro': _material_specs(cls = hydro,
                            type = 7,
                            nature = 7)
}
