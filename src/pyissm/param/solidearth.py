from .. import utils
from . import param_utils
from . import class_registry
from . import rotational
from . import lovenumbers
from .. import execute

## ------------------------------------------------------
## solidearth.earth
## ------------------------------------------------------
@class_registry.register_class
class earth(class_registry.manage_state):
    """
    Earth solid earth configuration class for ISSM.

    This class configures the solid earth response model for Earth, including
    glacial isostatic adjustment (GIA), sea level changes, and rotational feedbacks.
    It handles the coupling between ice sheet dynamics and solid earth deformation
    through Love numbers, rotational parameters, and partition vectors.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    settings : settings
        Solid earth settings object containing solver parameters and physical options.
    external : any, default=None
        External solution of the type solidearthsolution.
    lovenumbers : lovenumbers
        Love numbers object defining elastic and viscous earth response.
    rotational : rotational
        Rotational parameters object for polar motion calculations.
    planetradius : float
        Earth's radius [m]. Automatically set using utils.planetradius('earth').
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested from the solid earth model.
    transfercount : str, default='List of transfer count'
        Number of ice caps that vertices are part of.
    transitions : str, default='List of transitions'
        Indices into parts of the mesh that will be ice caps.
    partitionice : str, default='List of partionice'
        Ice partition vector for barystatic contribution.
    partitionhydro : str, default='List of partitionhydro'
        Hydro partition vector for barystatic contribution.
    partitionocean : str, default='List of partitionocean'
        Ocean partition vector for barystatic contribution.

    Methods
    -------
    __init__(self, other=None)
        Initializes the Earth solid earth parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the solid earth parameters.
    __str__(self)
        Returns a short string identifying the class.

    Notes
    -----
    This class is specifically configured for Earth with the appropriate planetary
    radius. For other planetary bodies, use the corresponding classes (e.g., europa).
    
    The solid earth model computes glacial isostatic adjustment effects including:
    - Elastic and viscous deformation
    - Self-attraction and loading effects
    - Rotational feedbacks
    - Sea level redistribution

    Examples
    --------
    md.solidearth = pyissm.param.solidearth.earth()
    md.solidearth.settings.elastic = 1
    md.solidearth.settings.viscous = 1
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.settings          = settings()
        self.external          = None
        self.lovenumbers       = lovenumbers()
        self.rotational        = rotational()
        self.planetradius      = utils.planetradius('earth')
        self.requested_outputs = 'List of requested outputs'
        self.transfercount     = 'List of transfer count'
        self.transitions       = 'List of transitions'
        self.partitionice      = 'List of partionice'
        self.partitionhydro    = 'List of partitionhydro'
        self.partitionocean    = 'List of partitionocean'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   solidearth inputs, forcings, and settings:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'planetradius', 'planet radius [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'transitions', 'indices into parts of the mesh that will be icecaps'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'transfercount', 'number of icecaps vertices are part of'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'partitionice', 'ice partition vector for barystatic contribution'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'partitionhydro', 'hydro partition vector for barystatic contribution'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'partitionocean', 'ocean partition vector for barystatic contribution'))
        if not self.external:
            s += '{}\n'.format(param_utils.fielddisplay(self, 'external', 'external solution, of the type solidearthsolution'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - solidearth.earth Class'
        return s

    # Marshall method for saving the earth parameters
    def marshall_class(self, prefix, md, fid):
        """
        Marshall the earth parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.

        Returns
        -------
        None
        """
        self.settings.marshall_class(prefix, md, fid)

        ## TODO: Marshall requested outputs and add all other fields


## ------------------------------------------------------
## solidearth.europa
## ------------------------------------------------------
@class_registry.register_class
class europa(class_registry.manage_state):
    """
    Europa solid earth configuration class for ISSM.

    This class configures the solid earth response model for Jupiter's moon Europa,
    including tidal heating, ice shell dynamics, and subsurface ocean interactions.
    It handles the coupling between ice dynamics and solid body deformation with
    parameters appropriate for Europa's smaller size and different composition.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    settings : settings
        Solid earth settings object containing solver parameters and physical options.
    external : any, default=None
        External solution of the type solidearthsolution.
    lovenumbers : lovenumbers
        Love numbers object defining elastic and viscous body response for Europa.
    rotational : rotational
        Rotational parameters object for Europa's rotation and tidal effects.
    planetradius : float
        Europa's radius [m]. Automatically set using utils.planetradius('europa').
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested from the solid body model.
    transfercount : str, default='List of transfer count'
        Number of ice caps that vertices are part of.
    transitions : str, default='List of transitions'
        Indices into parts of the mesh that will be ice caps.
    partitionice : str, default='List of partionice'
        Ice partition vector for barystatic contribution.
    partitionhydro : str, default='List of partitionhydro'
        Hydro partition vector for barystatic contribution.
    partitionocean : str, default='List of partitionocean'
        Ocean partition vector for barystatic contribution.

    Methods
    -------
    __init__(self, other=None)
        Initializes the Europa solid body parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the solid body parameters.
    __str__(self)
        Returns a short string identifying the class.

    Notes
    -----
    This class is specifically configured for Europa with the appropriate planetary
    radius and physical parameters. Europa's smaller size and different composition
    require different Love numbers and rotational parameters compared to Earth.
    
    The solid body model for Europa accounts for:
    - Tidal heating and deformation
    - Ice shell-ocean coupling
    - Smaller gravitational field effects
    - Different material properties

    Examples
    --------
    md.solidearth = pyissm.param.solidearth.europa()
    md.solidearth.settings.elastic = 1
    md.solidearth.lovenumbers = custom_europa_lovenumbers
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.settings          = settings()
        self.external          = None
        self.lovenumbers       = lovenumbers()
        self.rotational        = rotational()
        self.planetradius      = utils.planetradius('europa')
        self.requested_outputs = 'List of requested outputs'
        self.transfercount     = 'List of transfer count'
        self.transitions       = 'List of transitions'
        self.partitionice      = 'List of partionice'
        self.partitionhydro    = 'List of partitionhydro'
        self.partitionocean    = 'List of partitionocean'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   solidearth inputs, forcings, and settings:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'planetradius', 'planet radius [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'transitions', 'indices into parts of the mesh that will be icecaps'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'transfercount', 'number of icecaps vertices are part of'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'partitionice', 'ice partition vector for barystatic contribution'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'partitionhydro', 'hydro partition vector for barystatic contribution'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'partitionocean', 'ocean partition vector for barystatic contribution'))
        if not self.external:
            s += '{}\n'.format(param_utils.fielddisplay(self, 'external', 'external solution, of the type solidearthsolution'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - solidearth.europa Class'
        return s

## ------------------------------------------------------
## solidearth.settings
## ------------------------------------------------------
@class_registry.register_class
class settings(class_registry.manage_state):
    """
    Solid earth solver settings and physical options for ISSM.

    This class contains all the configuration parameters for the solid earth solver,
    including convergence criteria, physical process toggles, numerical accuracy
    settings, and model type selection for glacial isostatic adjustment calculations.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    reltol : float, default=0
        Sea level change relative convergence criterion. Default 0: not applied.
    abstol : float, default=0
        Sea level change absolute convergence criterion. Default 0: not applied.
    maxiter : int, default=0
        Maximum number of nonlinear iterations.
    selfattraction : int, default=1
        Enables surface mass load to perturb the gravity field (0 or 1).
    elastic : int, default=1
        Enables elastic deformation from surface loading (0 or 1).
    viscous : int, default=1
        Enables viscous deformation from surface loading (0 or 1).
    rotation : int, default=1
        Enables polar motion to feedback on the GRD fields (0 or 1).
    grdocean : int, default=1
        Does this planet have an ocean. If 1: global water mass is conserved in GRD module.
    ocean_area_scaling : float, default=0
        Correction for model representation of ocean area. Default 0: no correction.
    runfrequency : int, default=1
        How many time steps to skip before running solid earth solver during transient.
    sealevelloading : int, default=1
        Enables surface loading from sea-level change (0 or 1).
    isgrd : int, default=0
        Compute GRD patterns (0 or 1). Default 0: not computed.
    compute_bp_grd : int, default=0
        Compute GRD patterns for bottom pressure loads (0 or 1). Default 0: not computed.
    degacc : float, default=0
        Accuracy for numerical discretization of Green's functions [degrees]. Default: 0.01 deg.
    timeacc : float, default=0
        Time accuracy for numerical discretization of Green's functions [years]. Default: 1 year.
    horiz : int, default=0
        Horizontal displacement calculation flag.
    grdmodel : int, default=0
        Type of deformation model. 0: no GRD, 1: spherical GRD (SESAW), 2: half-space planar GRD (Ivins).
    cross_section_shape : int, default=0
        Cross-section shape for loading. 1: square-edged (default), 2: elliptical.

    Methods
    -------
    __init__(self, other=None)
        Initializes the solid earth settings, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the solid earth settings.
    __str__(self)
        Returns a short string identifying the class.

    Notes
    -----
    The solid earth solver supports multiple deformation models:
    - grdmodel=0: No glacial isostatic adjustment
    - grdmodel=1: Spherical earth model (SESAW approach)
    - grdmodel=2: Half-space planar model (visco-elastic, Ivins approach)

    The solver can be configured to include various physical processes:
    - Elastic deformation (instantaneous response)
    - Viscous deformation (time-dependent response)  
    - Self-attraction and loading effects
    - Rotational feedbacks and polar motion

    Examples
    --------
    md.solidearth.settings = pyissm.param.solidearth.settings()
    md.solidearth.settings.elastic = 1
    md.solidearth.settings.viscous = 1
    md.solidearth.settings.grdmodel = 1
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.reltol = 0
        self.abstol = 0
        self.maxiter = 0
        self.selfattraction = 1
        self.elastic = 1
        self.viscous = 1
        self.rotation = 1
        self.grdocean = 1
        self.ocean_area_scaling = 0
        self.runfrequency = 1
        self.sealevelloading = 1
        self.isgrd = 0
        self.compute_bp_grd = 0
        self.degacc = 0
        self.timeacc = 0
        self.horiz = 0
        self.grdmodel = 0
        self.cross_section_shape = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   solidearth settings:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'reltol', 'sea level change relative convergence criterion (default, NaN: not applied)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'abstol', 'sea level change absolute convergence criterion (default, NaN: not applied)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'maxiter', 'maximum number of nonlinear iterations'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'grdocean', 'does this planet have an ocean, if set to 1: global water mass is conserved in GRD module (default: 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ocean_area_scaling', 'correction for model representation of ocean area (default: No correction)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sealevelloading', 'enables surface loading from sea-level change (default: 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isgrd', 'compute GRD patterns (default: 1'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'compute_bp_grd', 'compute GRD patterns for bottom pressure loads (default 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'runfrequency', 'how many time steps we skip before we run solidearthsettings solver during transient (default: 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'selfattraction', 'enables surface mass load to perturb the gravity field'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'elastic', 'enables elastic deformation from surface loading'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'viscous', 'enables viscous deformation from surface loading'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rotation', 'enables polar motion to feedback on the GRD fields'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'degacc', 'accuracy (default: .01 deg) for numerical discretization of the Green\'s functions'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'timeacc', 'time accuracy (default: 1 year) for numerical discretization of the Green\'s functions'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'grdmodel', 'type of deformation model, 0 for no GRD, 1 for spherical GRD model (SESAW model), 2 for half-space planar GRD (visco-elastic model from Ivins)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'cross_section_shape', '1: square-edged (default). 2: elliptical. See iedge in GiaDeflectionCore'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - solidearth.settings Class'
        return s

    
   # Marshall method for saving the settings parameters
    def marshall_class(self, prefix, md, fid):
        """
        Marshall the settings parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.

        Returns
        -------
        None
        """
        ## Write Double fields
        fieldnames = ['reltol', 'abstol', 'degacc']
        for fieldname in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = fieldname, format = 'Double', name = 'md.solidearth.settings.' + fieldname)

        ## Write Integer fields
        fieldnames = ['maxiter', 'runfrequency', 'horiz', 'sealevelloading', 'isgrd', 'compute_bp_grd', 'grdmodel', 'cross_section_shape']
        for fieldname in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = fieldname, format = 'Integer', name = 'md.solidearth.settings.' + fieldname)
        ## Write Boolean fields
        fieldnames = ['selfattraction', 'elastic', 'viscous', 'rotation', 'grdocean', 'ocean_area_scaling']
        for fieldname in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = fieldname, format = 'Boolean', name = 'md.solidearth.settings.' + fieldname)
        
        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'timeacc', format = 'Double', name = 'md.solidearth.settings.timeacc', scale = md.constants.yts)


## ------------------------------------------------------
## solidearth.solution
## ------------------------------------------------------
@class_registry.register_class
class solution(class_registry.manage_state):
    """
    Solid earth solution container class for ISSM.

    This class stores the time series results from solid earth deformation
    calculations, including bedrock displacement components and geoid changes.
    All solutions are time series data with units specified for each field.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    displacementeast : str, default='displacementeast timeseries'
        Solid-earth eastwards bedrock displacement time series [m].
    displacementnorth : str, default='displacementsnorth timeseries'  
        Solid-earth northwards bedrock displacement time series [m].
    displacementup : str, default='displacementup timeseries'
        Solid-earth bedrock uplift time series [m].
    geoid : str, default='geoid timeseries'
        Solid-earth geoid time series [m].

    Methods
    -------
    __init__(self, other=None)
        Initializes the solution container, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the solution fields.
    __str__(self)
        Returns a short string identifying the class.

    Notes
    -----
    All time series have units of years (yr) for time and meters (m) for displacements
    and geoid heights. The solution represents the solid earth response to surface
    loading from ice mass changes.
    
    The displacement fields represent:
    - displacementeast: Horizontal movement in eastward direction
    - displacementnorth: Horizontal movement in northward direction  
    - displacementup: Vertical movement (positive upward)
    - geoid: Changes in geoid height due to mass redistribution

    Examples
    --------
    # Solution is typically populated by the solver
    sol = md.solidearth.solution
    east_disp = sol.displacementeast
    uplift = sol.displacementup
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.displacementeast = 'displacementeast timeseries'
        self.displacementnorth = 'displacementsnorth timeseries'
        self.displacementup = 'displacementup timeseries'
        self.geoid = 'geoid timeseries'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   solidearth solution:\n'
        s += '         units for time series is (yr)\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'displacementeast', 'solid-Earth Eastwards bedrock displacement series (m)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'displacementnorth', 'solid-Earth Northwards bedrock displacement time series (m)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'displacementup', 'solid-Earth bedrock uplift time series (m)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'geoid', 'solid-Earth geoid time series (m)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - solidearth.solution Class'
        return s