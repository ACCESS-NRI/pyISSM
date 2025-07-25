import numpy as np
from . import param_utils
from . import class_registry

## ------------------------------------------------------
## smb.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    """
    Default surface mass balance (SMB) parameters class for ISSM.

    This class encapsulates the default parameters for surface mass balance in the ISSM (Ice Sheet System Model) framework.
    It defines the main SMB-related parameters.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    mass_balance : ndarray, default=np.nan
        Surface mass balance [m/yr ice eq].
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested (default: ['SmbMassBalance']).
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic (default), 1: Geometric, 2: Harmonic.

    Methods
    -------
    __init__(self, other=None)
        Initializes the SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the SMB parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.smb = pyissm.param.smb.default()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.mass_balance = np.nan
        self.steps_per_step = 1
        self.requested_outputs = 'List of requested outputs' # Default = ['SmbMassBalance']
        self.averaging = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'mass_balance', 'surface mass balance [m/yr ice eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.default Class'
        return s

## ------------------------------------------------------
## smb.arma
## ------------------------------------------------------
@class_registry.register_class
class arma(class_registry.manage_state):
    """
    ARMA (AutoRegressive Moving Average) surface mass balance model for ISSM.

    This class implements an ARMA-based surface mass balance model that combines 
    autoregressive and moving average components with piecewise polynomial trends
    and elevation-dependent lapse rates for basin-specific SMB modeling.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    num_basins : int, default=0
        Number of different basins [unitless].
    num_params : int, default=0
        Number of different parameters in the piecewise-polynomial (1:intercept only, 
        2:with linear trend, 3:with quadratic trend, etc.).
    num_breaks : int, default=0
        Number of different breakpoints in the piecewise-polynomial (separating 
        num_breaks+1 periods).
    polynomialparams : ndarray, default=np.nan
        Coefficients for the polynomial (const,trend,quadratic,etc.), dim1 for basins,
        dim2 for periods, dim3 for orders.
    ar_order : float, default=0.0
        Order of the autoregressive model [unitless].
    ma_order : float, default=0.0
        Order of the moving-average model [unitless].
    arlag_coefs : ndarray, default=np.nan
        Basin-specific vectors of AR lag coefficients [unitless].
    malag_coefs : ndarray, default=np.nan
        Basin-specific vectors of MA lag coefficients [unitless].
    datebreaks : ndarray, default=np.nan
        Dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr].
    basin_id : ndarray, default=np.nan
        Basin number assigned to each element [unitless].
    lapserates : ndarray, default=np.nan
        Basin-specific SMB lapse rates applied in each elevation bin, 1 row per basin,
        1 column per bin, dimension 3 can be of size 12 to prescribe monthly varying 
        values [m ice eq yr^-1 m^-1].
    elevationbins : ndarray, default=np.nan
        Basin-specific separations between elevation bins, 1 row per basin, 1 column 
        per limit between bins [m].
    refelevation : ndarray, default=np.nan
        Basin-specific reference elevations at which SMB is calculated [m].
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the ARMA SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the ARMA SMB parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.smb = pyissm.param.smb.arma()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.num_basins = 0
        self.num_params = 0
        self.num_breaks = 0
        self.polynomialparams = np.nan
        self.ar_order = 0.0
        self.ma_order = 0.0
        self.arlag_coefs = np.nan
        self.malag_coefs = np.nan
        self.polynomialparams = np.nan
        self.datebreaks = np.nan
        self.basin_id = np.nan
        self.lapserates = np.nan
        self.elevationbins = np.nan
        self.refelevation = np.nan
        self.datebreaks = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'num_basins', 'number of different basins [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'basin_id', 'basin number assigned to each element [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'num_breaks', 'number of different breakpoints in the piecewise-polynomial (separating num_breaks+1 periods)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'num_params', 'number of different parameters in the piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'polynomialparams', 'coefficients for the polynomial (const,trend,quadratic,etc.),dim1 for basins,dim2 for periods,dim3 for orders, ex: polyparams=cat(num_params,intercepts,trendlinearcoefs,trendquadraticcoefs)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'datebreaks', 'dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ar_order', 'order of the autoregressive model [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ma_order', 'order of the moving-average model [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'arma_timestep', 'time resolution of the ARMA model [yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'arlag_coefs', 'basin-specific vectors of AR lag coefficients [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'malag_coefs', 'basin-specific vectors of MA lag coefficients [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'lapserates', 'basin-specific SMB lapse rates applied in each elevation bin, 1 row per basin, 1 column per bin, dimension 3 can be of size 12 to prescribe monthly varying values [m ice eq yr^-1 m^-1] (default: no lapse rate)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'elevationbins', 'basin-specific separations between elevation bins, 1 row per basin, 1 column per limit between bins, dimension 3 can be of size 12 to prescribe monthly varying values [m] (default: no basin separation)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'refelevation', 'basin-specific reference elevations at which SMB is calculated, and from which SMB is downscaled using lapserates (default: basin mean elevation) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.arma Class'
        return s

## ------------------------------------------------------
## smb.components
## ------------------------------------------------------
@class_registry.register_class
class components(class_registry.manage_state):
    """
    Component-based surface mass balance model for ISSM.

    This class implements a component-based SMB model where the surface mass balance
    is calculated as SMB = accumulation - runoff - evaporation. Each component can
    be specified independently.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    accumulation : ndarray, default=np.nan
        Accumulated snow [m/yr ice eq].
    runoff : ndarray, default=np.nan
        Amount of ice melt lost from the ice column [m/yr ice eq].
    evaporation : ndarray, default=np.nan
        Amount of ice lost to evaporative processes [m/yr ice eq].
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the component SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the component SMB parameters.
    __str__(self)
        Returns a short string identifying the class.

    Notes
    -----
    The surface mass balance is computed as:
    SMB = accumulation - runoff - evaporation

    Examples
    --------
    md.smb = pyissm.param.smb.components()
    md.smb.accumulation = accumulation_data
    md.smb.runoff = runoff_data
    md.smb.evaporation = evaporation_data
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.accumulation = np.nan
        self.runoff = np.nan
        self.evaporation = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters (SMB=accumulation-runoff-evaporation) :\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'accumulation', 'accumulated snow [m/yr ice eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'runoff', 'amount of ice melt lost from the ice column [m/yr ice eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'evaporation', 'mount of ice lost to evaporative processes [m/yr ice eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.components Class'
        return s

## ------------------------------------------------------
## smb.d18opdd
## ------------------------------------------------------
@class_registry.register_class
class d18opdd(class_registry.manage_state):
    """
    Delta-18-O driven positive degree day surface mass balance model for ISSM.

    This class implements a positive degree day (PDD) SMB model driven by delta-18-O 
    isotope data for paleoclimate applications. It includes temperature and precipitation
    scaling based on isotope ratios and elevation-dependent corrections.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    desfac : float, default=0.5
        Desertification elevation factor (between 0 and 1) [m].
    s0p : ndarray, default=np.nan
        Elevation from precipitation source (between 0 and a few 1000s m) [m].
    s0t : ndarray, default=np.nan
        Elevation from temperature source (between 0 and a few 1000s m) [m].
    rlaps : float, default=6.5
        Present day lapse rate [degree/km].
    rlapslgm : float, default=6.5
        LGM lapse rate [degree/km].
    dpermil : float, default=2.4
        Degree per mil, required if d18opd is activated.
    f : float, default=0.169
        Precipitation/temperature scaling factor, required if d18opd is activated.
    Tdiff : ndarray, default=np.nan
        Temperature difference field.
    sealev : ndarray, default=np.nan
        Sea level data.
    ismungsm : int, default=0
        Is mungsm parametrisation activated (0 or 1).
    isd18opd : int, default=1
        Is delta18o parametrisation from present day temperature and precipitation activated (0 or 1).
    issetpddfac : int, default=0
        Is user passing in defined PDD factors (0 or 1).
    istemperaturescaled : int, default=1
        Is temperature scaled to delta18o value (0 or 1).
    isprecipscaled : int, default=1
        Is precipitation scaled to delta18o value (0 or 1).
    delta18o : ndarray, default=np.nan
        Delta-18-O values [per mil].
    delta18o_surface : ndarray, default=np.nan
        Surface delta-18-O values.
    temperatures_presentday : ndarray, default=np.nan
        Monthly present day surface temperatures [K].
    precipitations_presentday : ndarray, default=np.nan
        Monthly surface precipitation [m/yr water eq].
    temperatures_reconstructed : ndarray, default=np.nan
        Monthly historical surface temperatures [K].
    precipitations_reconstructed : ndarray, default=np.nan
        Monthly historical precipitation [m/yr water eq].
    pddfac_snow : ndarray, default=np.nan
        PDD factor for snow [mm ice equiv/day/degree C].
    pddfac_ice : ndarray, default=np.nan
        PDD factor for ice [mm ice equiv/day/degree C].
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the delta-18-O PDD SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the delta-18-O PDD SMB parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.smb = pyissm.param.smb.d18opdd()
    md.smb.delta18o = delta18o_data
    md.smb.temperatures_presentday = temp_data
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.desfac = 0.5
        self.s0p = np.nan
        self.s0t = np.nan
        self.rlaps = 6.5
        self.rlapslgm = 6.5
        self.dpermil = 2.4
        self.f = 0.169
        self.Tdiff = np.nan
        self.sealev = np.nan
        self.ismungsm = 0
        self.isd18opd = 1
        self.issetpddfac = 0
        self.istemperaturescaled = 1
        self.isprecipscaled = 1
        self.delta18o = np.nan
        self.delta18o_surface = np.nan
        self.temperatures_presentday = np.nan
        self.precipitations_presentday = np.nan
        self.temperatures_reconstructed = np.nan
        self.precipitations_reconstructed = np.nan
        self.pddfac_snow = np.nan
        self.pddfac_ice = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'isd18opd', 'is delta18o parametrisation from present day temperature and precipitation activated (0 or 1, default is 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'issetpddfac', 'is user passing in defined pdd factors (0 or 1, default is 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'desfac', 'desertification elevation factor (between 0 and 1, default is 0.5) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 's0p', 'should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 's0t', 'should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rlaps', 'present day lapse rate [degree/km]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'temperatures_presentday', 'monthly present day surface temperatures [K], required if delta18o/mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'precipitations_presentday', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'istemperaturescaled', 'if delta18o parametrisation from present day temperature and precipitation is activated, is temperature scaled to delta18o value? (0 or 1, default is 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isprecipscaled', 'if delta18o parametrisation from present day temperature and precipitation is activated, is precipitation scaled to delta18o value? (0 or 1, default is 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'temperatures_reconstructed', 'monthly historical surface temperatures [K], required if delta18o/mungsm/d18opd is activated and istemperaturescaled is not activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'precipitations_reconstructed', 'monthly historical precipitation [m/yr water eq], required if delta18o/mungsm/d18opd is activated and isprecipscaled is not activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'delta18o', 'delta18o [per mil], required if pdd is activated and delta18o activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'dpermil', 'degree per mil, required if d18opd is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'f', 'precip/temperature scaling factor, required if d18opd is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'pddfac_snow', 'Pdd factor for snow for all the domain [mm ice equiv/day/degree C]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'pddfac_ice', 'Pdd factor for ice for all the domain [mm ice equiv/day/degree C]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.d18opdd Class'
        return s

## ------------------------------------------------------
## smb.gradients
## ------------------------------------------------------
@class_registry.register_class
class gradients(class_registry.manage_state):
    """
    Gradient-based surface mass balance model for ISSM.

    This class implements a gradient-based SMB model where SMB varies linearly with
    elevation relative to a reference elevation and SMB. Different gradients can be
    specified for accumulation and ablation regimes.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    href : ndarray, default=np.nan
        Reference elevation from which deviation is used to calculate SMB adjustment [m].
    smbref : ndarray, default=np.nan
        Reference SMB from which deviation is calculated [m/yr ice equiv].
    b_pos : ndarray, default=np.nan
        Slope of elevation-SMB regression line for accumulation regime.
    b_neg : ndarray, default=np.nan
        Slope of elevation-SMB regression line for ablation regime.
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the gradient SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the gradient SMB parameters.
    __str__(self)
        Returns a short string identifying the class.

    Notes
    -----
    SMB is calculated as:
    SMB = smbref + gradient * (elevation - href)
    where gradient = b_pos for positive SMB or b_neg for negative SMB.

    Examples
    --------
    md.smb = pyissm.param.smb.gradients()
    md.smb.href = reference_elevation
    md.smb.smbref = reference_smb
    md.smb.b_pos = positive_gradient
    md.smb.b_neg = negative_gradient
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.href = np.nan
        self.smbref = np.nan
        self.b_pos = np.nan
        self.b_neg = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'issmbgradients', 'is smb gradients method activated (0 or 1, default is 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'href', 'reference elevation from which deviation is used to calculate SMB adjustment in smb gradients method'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'smbref', 'reference smb from which deviation is calculated in smb gradients method [m/yr ice equiv]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'b_pos', 'slope of hs - smb regression line for accumulation regime required if smb gradients is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'b_neg', 'slope of hs - smb regression line for ablation regime required if smb gradients is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.gradients Class'
        return s

## ------------------------------------------------------
## smb.gradientscomponents
## ------------------------------------------------------
@class_registry.register_class
class gradientscomponents(class_registry.manage_state):
    """
    Component-based gradient surface mass balance model for ISSM.

    This class implements a gradient-based SMB model where accumulation and runoff
    components vary separately with elevation. Each component has its own reference
    value, reference elevation, and elevation gradient.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    accuref : ndarray, default=np.nan
        Reference value of the accumulation [m ice eq/yr].
    accualti : ndarray, default=np.nan
        Altitude at which the accumulation is equal to the reference value [m].
    accugrad : ndarray, default=np.nan
        Gradient of the variation of the accumulation (0 for uniform accumulation) [m ice eq/yr/m].
    runoffref : ndarray, default=np.nan
        Reference value of the runoff [m w.e. y^-1].
    runoffalti : ndarray, default=np.nan
        Altitude at which the runoff is equal to the reference value [m].
    runoffgrad : ndarray, default=np.nan
        Gradient of the variation of the runoff (0 for uniform runoff) [m w.e. m^-1 y^-1].
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the gradient components SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the gradient components SMB parameters.
    __str__(self)
        Returns a short string identifying the class.

    Notes
    -----
    SMB components are calculated as:
    accumulation = accuref + accugrad * (elevation - accualti)
    runoff = runoffref + runoffgrad * (elevation - runoffalti)
    SMB = accumulation - runoff

    Examples
    --------
    md.smb = pyissm.param.smb.gradientscomponents()
    md.smb.accuref = reference_accumulation
    md.smb.runoffref = reference_runoff
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.accuref = np.nan
        self.accualti = np.nan
        self.accugrad = np.nan
        self.runoffref = np.nan
        self.runoffalti = np.nan
        self.runoffgrad = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'issmbgradients', 'is smb gradients method activated (0 or 1, default is 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'accuref', ' reference value of the accumulation'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'accualti', ' Altitude at which the accumulation is equal to the reference value'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'accugrad', ' Gradient of the variation of the accumulation (0 for uniform accumulation)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'runoffref', ' reference value of the runoff m w.e. y-1 (temperature times ddf)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'runoffalti', ' Altitude at which the runoff is equal to the reference value'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'runoffgrad', ' Gradient of the variation of the runoff (0 for uniform runoff) m w.e. m-1 y-1 (lapse rate times ddf)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.gradientscomponents Class'
        return s

## ------------------------------------------------------
## smb.gradientsela
## ------------------------------------------------------
@class_registry.register_class
class gradientsela(class_registry.manage_state):
    """
    Equilibrium Line Altitude (ELA) gradient surface mass balance model for ISSM.

    This class implements an ELA-based SMB model where SMB varies linearly with 
    elevation relative to the equilibrium line altitude. Different gradients are
    applied above and below the ELA, with optional caps on maximum and minimum SMB rates.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    ela : ndarray, default=np.nan
        Equilibrium line altitude from which deviation is used to calculate SMB [m a.s.l.].
    b_pos : ndarray, default=np.nan
        Vertical SMB gradient (dB/dz) above ELA [m ice eq./yr/m].
    b_neg : ndarray, default=np.nan
        Vertical SMB gradient (dB/dz) below ELA [m ice eq./yr/m].
    b_max : float, default=9999
        Upper cap on SMB rate [m ice eq./yr]. Default: 9999 (no cap).
    b_min : float, default=-9999
        Lower cap on SMB rate [m ice eq./yr]. Default: -9999 (no cap).
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the ELA gradient SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the ELA gradient SMB parameters.
    __str__(self)
        Returns a short string identifying the class.

    Notes
    -----
    SMB is calculated as:
    - For elevation > ELA: SMB = b_pos * (elevation - ELA)
    - For elevation < ELA: SMB = b_neg * (elevation - ELA)
    SMB is then clamped between b_min and b_max if specified.

    Examples
    --------
    md.smb = pyissm.param.smb.gradientsela()
    md.smb.ela = equilibrium_line_altitude
    md.smb.b_pos = positive_gradient  # Above ELA
    md.smb.b_neg = negative_gradient  # Below ELA
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.ela = np.nan
        self.b_pos = np.nan
        self.b_neg = np.nan
        self.b_max = 9999
        self.b_min = -9999
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

        s += '\n   SMB gradients ela parameters:'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ela', ' equilibrium line altitude from which deviation is used to calculate smb using the smb gradients ela method [m a.s.l.]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'b_pos', ' vertical smb gradient (dB/dz) above ela'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'b_neg', ' vertical smb gradient (dB/dz) below ela'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'b_max', ' upper cap on smb rate, default: 9999 (no cap) [m ice eq./yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'b_min', ' lower cap on smb rate, default: -9999 (no cap) [m ice eq./yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.gradientsela Class'
        return s

## ------------------------------------------------------
## smb.henning
## ------------------------------------------------------
@class_registry.register_class
class henning(class_registry.manage_state):
    """
    Henning surface mass balance model for ISSM.

    This class implements the Henning SMB parametrization, which is a specialized
    approach for modeling surface mass balance in ice sheet applications.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    smbref : ndarray, default=np.nan
        Reference surface mass balance [m/yr ice eq].
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the Henning SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the Henning SMB parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.smb = pyissm.param.smb.henning()
    md.smb.smbref = reference_smb_data
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.smbref = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'mass_balance', 'surface mass balance [m/yr ice eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.henning Class'
        return s

## ------------------------------------------------------
## smb.meltcomponents
## ------------------------------------------------------
@class_registry.register_class
class meltcomponents(class_registry.manage_state):
    """
    Melt component-based surface mass balance model for ISSM.

    This class implements a component-based SMB model that explicitly separates
    melt and refreeze processes. The surface mass balance is calculated as 
    SMB = accumulation - evaporation - melt + refreeze.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    accumulation : ndarray, default=np.nan
        Accumulated snow [m/yr ice eq].
    evaporation : ndarray, default=np.nan
        Amount of ice lost to evaporative processes [m/yr ice eq].
    melt : ndarray, default=np.nan
        Amount of ice melt in the ice column [m/yr ice eq].
    refreeze : ndarray, default=np.nan
        Amount of ice melt refrozen in the ice column [m/yr ice eq].
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the melt components SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the melt components SMB parameters.
    __str__(self)
        Returns a short string identifying the class.

    Notes
    -----
    The surface mass balance is computed as:
    SMB = accumulation - evaporation - melt + refreeze

    This formulation explicitly accounts for refreezing processes that can occur
    in firn layers, which is important for accurate SMB modeling in cold regions.

    Examples
    --------
    md.smb = pyissm.param.smb.meltcomponents()
    md.smb.accumulation = accumulation_data
    md.smb.melt = melt_data
    md.smb.refreeze = refreeze_data
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.accumulation = np.nan
        self.evaporation = np.nan
        self.melt = np.nan
        self.refreeze = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters with melt (SMB = accumulation-evaporation-melt+refreeze):\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'accumulation', 'accumulated snow [m/yr ice eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'evaporation', 'mount of ice lost to evaporative processes [m/yr ice eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'melt', 'amount of ice melt in the ice column [m/yr ice eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'refreeze', 'amount of ice melt refrozen in the ice column [m/yr ice eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.meltcomponents Class'
        return s

## ------------------------------------------------------
## smb.pdd
## ------------------------------------------------------
@class_registry.register_class
class pdd(class_registry.manage_state):
    """
    Positive Degree Day surface mass balance model for ISSM.

    This class implements a positive degree day (PDD) SMB model that calculates
    surface mass balance based on temperature and precipitation data. It supports
    multiple temperature and precipitation data sources, including delta-18-O and
    MUNGSM parametrizations for paleoclimate applications.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    precipitation : ndarray, default=np.nan
        Monthly surface precipitation [m/yr water eq].
    monthlytemperatures : ndarray, default=np.nan
        Monthly surface temperatures [K].
    desfac : float, default=0.5
        Desertification elevation factor (between 0 and 1) [m].
    s0p : ndarray, default=np.nan
        Elevation from precipitation source (between 0 and a few 1000s m) [m].
    s0t : ndarray, default=np.nan
        Elevation from temperature source (between 0 and a few 1000s m) [m].
    rlaps : float, default=6.5
        Present day lapse rate [degree/km].
    rlapslgm : float, default=6.5
        LGM lapse rate [degree/km].
    Pfac : ndarray, default=np.nan
        Time interpolation parameter for precipitation, 1D(year).
    Tdiff : ndarray, default=np.nan
        Time interpolation parameter for temperature, 1D(year).
    sealev : ndarray, default=np.nan
        Sea level [m], 1D(year).
    isdelta18o : int, default=0
        Is temperature and precipitation delta18o parametrisation activated (0 or 1).
    ismungsm : int, default=0
        Is temperature and precipitation mungsm parametrisation activated (0 or 1).
    issetpddfac : int, default=0
        Is user passing in defined PDD factors (0 or 1).
    delta18o : float, default=0
        Delta-18-O values [per mil].
    delta18o_surface : ndarray, default=np.nan
        Surface elevation of the delta18o site [m].
    temperatures_presentday : ndarray, default=np.nan
        Monthly present day surface temperatures [K].
    temperatures_lgm : ndarray, default=np.nan
        Monthly LGM surface temperatures [K].
    precipitations_presentday : ndarray, default=np.nan
        Monthly present day surface precipitation [m/yr water eq].
    precipitations_lgm : ndarray, default=np.nan
        Monthly LGM surface precipitation [m/yr water eq].
    pddfac_snow : ndarray, default=np.nan
        PDD factor for snow [mm ice equiv/day/degree C].
    pddfac_ice : ndarray, default=np.nan
        PDD factor for ice [mm ice equiv/day/degree C].
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the PDD SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the PDD SMB parameters.
    __str__(self)
        Returns a short string identifying the class.

    Notes
    -----
    The PDD model calculates melt based on the number of positive degree days,
    which is the sum of temperatures above freezing over a given time period.
    This approach is widely used in glaciology for its simplicity and effectiveness.

    Examples
    --------
    md.smb = pyissm.param.smb.pdd()
    md.smb.monthlytemperatures = temperature_data
    md.smb.precipitation = precipitation_data
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.precipitation = np.nan
        self.monthlytemperatures = np.nan
        self.desfac = 0.5
        self.s0p = np.nan
        self.s0t = np.nan
        self.rlaps = 6.5
        self.rlapslgm = 6.5
        self.Pfac = np.nan
        self.Tdiff = np.nan
        self.sealev = np.nan
        self.isdelta18o = 0
        self.ismungsm = 0
        self.issetpddfac = 0
        self.delta18o = 0
        self.delta18o_surface = np.nan
        self.temperatures_presentday = np.nan
        self.temperatures_lgm = np.nan
        self.precipitations_presentday = np.nan
        self.precipitations_lgm = np.nan
        self.pddfac_snow = np.nan
        self.pddfac_ice = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'isdelta18o', 'is temperature and precipitation delta18o parametrisation activated (0 or 1, default is 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ismungsm', 'is temperature and precipitation mungsm parametrisation activated (0 or 1, default is 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'issetpddfac', 'is user passing in defined pdd factors (0 or 1, default is 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'desfac', 'desertification elevation factor (between 0 and 1, default is 0.5) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 's0p', 'should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 's0t', 'should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rlaps', 'present day lapse rate [degree/km]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rlapslgm', 'LGM lapse rate [degree/km]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'monthlytemperatures', 'monthly surface temperatures [K], required if pdd is activated and delta18o not activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'precipitation', 'monthly surface precipitation [m/yr water eq], required if pdd is activated and delta18o or mungsm not activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'delta18o', 'delta18o [per mil], required if pdd is activated and delta18o activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'delta18o_surface', 'surface elevation of the delta18o site, required if pdd is activated and delta18o activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'temperatures_presentday', 'monthly present day surface temperatures [K], required if delta18o/mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'temperatures_lgm', 'monthly LGM surface temperatures [K], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'precipitations_presentday', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'precipitations_lgm', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Tdiff', 'time interpolation parameter for temperature, 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sealev', 'sea level [m], 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'temperatures_presentday', 'monthly present day surface temperatures [K], required if delta18o/mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'temperatures_lgm', 'monthly LGM surface temperatures [K], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'precipitations_presentday', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'precipitations_lgm', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Pfac', 'time interpolation parameter for precipitation, 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Tdiff', 'time interpolation parameter for temperature, 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sealev', 'sea level [m], 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.pdd Class'
        return s

## ------------------------------------------------------
## smb.pddSicopolis
## ------------------------------------------------------
@class_registry.register_class
class pddSicopolis(class_registry.manage_state):
    """
    SICOPOLIS-style Positive Degree Day surface mass balance model for ISSM.

    This class implements the SICOPOLIS PDD scheme (Calov & Greve, 2005) for 
    surface mass balance calculations. It includes temperature and precipitation
    anomalies, firn warming effects, and desertification corrections based on
    the SICOPOLIS ice sheet model approach.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    precipitation : ndarray, default=np.nan
        Monthly surface precipitation [m/yr water eq].
    monthlytemperatures : ndarray, default=np.nan
        Monthly surface temperatures [K].
    temperature_anomaly : ndarray, default=np.nan
        Anomaly to monthly reference temperature (additive) [K].
    precipitation_anomaly : ndarray, default=np.nan
        Anomaly to monthly precipitation (multiplicative, e.g. q = q0*exp(0.070458*DeltaT)) [unitless].
    smb_corr : ndarray, default=np.nan
        Correction of SMB after PDD call [m/a].
    desfac : float, default=-np.log(2.0)/1000
        Desertification elevation factor. Default: -log(2.0)/1000.
    s0p : ndarray, default=np.nan
        Elevation from precipitation source (between 0 and a few 1000s m) [m].
    s0t : ndarray, default=np.nan
        Elevation from temperature source (between 0 and a few 1000s m) [m].
    rlaps : float, default=7.4
        Present day lapse rate [degree/km]. Default: 7.4.
    isfirnwarming : int, default=1
        Is firn warming (Reeh 1991) activated (0 or 1). Default: 1.
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested (TemperaturePDD, SmbAccumulation, SmbMelt).

    Methods
    -------
    __init__(self, other=None)
        Initializes the SICOPOLIS PDD SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the SICOPOLIS PDD SMB parameters.
    __str__(self)
        Returns a short string identifying the class.

    Notes
    -----
    This implementation follows the SICOPOLIS PDD scheme as described in:
    Calov, R., & Greve, R. (2005). A semi-analytical solution for the positive 
    degree-day model with stochastic temperature variations. Journal of Glaciology, 
    51(172), 173-175.

    The firn warming correction (Reeh, 1991) adjusts melt rates based on firn
    temperature, which is important for accurate SMB calculations in accumulation zones.

    Examples
    --------
    md.smb = pyissm.param.smb.pddSicopolis()
    md.smb.monthlytemperatures = temperature_data
    md.smb.precipitation = precipitation_data
    md.smb.temperature_anomaly = temp_anomaly
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.precipitation = np.nan
        self.monthlytemperatures = np.nan
        self.temperature_anomaly = np.nan
        self.precipitation_anomaly = np.nan
        self.smb_corr = np.nan
        self.desfac = -np.log(2.0) / 1000
        self.s0p = np.nan
        self.s0t = np.nan
        self.rlaps = 7.4
        self.isfirnwarming = 1
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'
        s += '   SICOPOLIS PDD scheme (Calov & Greve, 2005):\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'monthlytemperatures', 'monthly surface temperatures [K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'precipitation', 'monthly surface precipitation [m/yr water eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'temperature_anomaly', 'anomaly to monthly reference temperature (additive [K])'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'precipitation_anomaly', 'anomaly to monthly precipitation (multiplicative, e.g. q = q0*exp(0.070458*DeltaT) after Huybrechts (2002)) [no unit])'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'smb_corr', 'correction of smb after PDD call [m/a]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 's0p', 'should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 's0t', 'should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rlaps', 'present day lapse rate (default is 7.4 degree/km)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'desfac', 'desertification elevation factor (default is -log(2.0)/1000)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isfirnwarming', 'is firnwarming (Reeh 1991) activated (0 or 1, default is 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested (TemperaturePDD, SmbAccumulation, SmbMelt)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.pddSicopolis Class'
        return s