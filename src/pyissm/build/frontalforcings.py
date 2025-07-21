import numpy as np
from . import build_utils
from . import class_registry

## ------------------------------------------------------
## frontalforcings.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    """
    Default frontalforcings parameters class for ISSM.

    This class encapsulates the default parameters for frontal forcings in the ISSM (Ice Sheet System Model) framework.
    It defines the main frontal forcing-related parameters.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    meltingrate : ndarray, default=np.nan
        Melting rate at given location [m/a].
    ablationrate : ndarray, default=np.nan
        Frontal ablation rate at given location [m/a], it contains both calving and melting.

    Methods
    -------
    __init__(self, other=None)
        Initializes the frontalforcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the frontalforcings parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.frontalforcings = pyissm.build.frontalforcings.default()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.meltingrate = np.nan
        self.ablationrate = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Frontalforcings parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'meltingrate', 'melting rate at given location [m/a]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'ablationrate', 'frontal ablation rate at given location [m/a], it contains both calving and melting'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - frontalforcings.default Class'
        return s

## ------------------------------------------------------
## frontalforcings.rignot
## ------------------------------------------------------
@class_registry.register_class
class rignot(class_registry.manage_state):
    """
    Rignot frontalforcings parameters class for ISSM.

    This class encapsulates the parameters for frontal forcings based on the Rignot methodology in the ISSM (Ice Sheet System Model) framework.
    It defines the main frontal forcing-related parameters specific to the Rignot approach.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    basin_id : ndarray, default=np.nan
        Basin ID for elements.
    num_basins : int, default=0
        Number of basins.
    subglacial_discharge : ndarray, default=np.nan
        Sum of subglacial discharge for each basin [m/d].
    thermalforcing : ndarray, default=np.nan
        Thermal forcing [°C].

    Methods
    -------
    __init__(self, other=None)
        Initializes the Rignot frontalforcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the Rignot frontalforcings parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.frontalforcings = pyissm.build.frontalforcings.rignot()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.basin_id = np.nan
        self.num_basins = 0
        self.subglacial_discharge = np.nan
        self.thermalforcing = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Frontalforcings parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'basin_id', 'basin ID for elements'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'num_basins', 'number of basins'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'subglacial_discharge', 'sum of subglacial discharge for each basin [m/d]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thermalforcing', 'thermal forcing [∘C]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - frontalforcings.rignot Class'
        return s

## ------------------------------------------------------
## frontalforcings.rignotarma
## ------------------------------------------------------
@class_registry.register_class
class rignotarma(class_registry.manage_state):
    """
    RignotARMA frontalforcings parameters class for ISSM.

    This class encapsulates the parameters for frontal forcings based on the Rignot methodology with ARMA (AutoRegressive Moving Average) modeling in the ISSM (Ice Sheet System Model) framework.
    It defines the main frontal forcing-related parameters specific to the RignotARMA approach, including polynomial trends, breakpoints, ARMA coefficients, and subglacial discharge modeling.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    num_basins : int, default=0
        Number of different basins.
    num_params : int, default=0
        Number of different parameters in the piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.).
    num_breaks : int, default=0
        Number of different breakpoints in the piecewise-polynomial (separating num_breaks+1 periods).
    polynomialparams : ndarray, default=np.nan
        Coefficients for the polynomial (const, trend, quadratic, etc.), dim1 for basins, dim2 for periods, dim3 for orders.
    datebreaks : ndarray, default=np.nan
        Dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr].
    ar_order : int, default=0
        Order of the autoregressive model.
    ma_order : int, default=0
        Order of the moving-average model.
    arma_timestep : int, default=0
        Time resolution of the ARMA model [yr].
    arlag_coefs : ndarray, default=np.nan
        Basin-specific vectors of AR lag coefficients.
    malag_coefs : ndarray, default=np.nan
        Basin-specific vectors of MA lag coefficients.
    monthlyvals_intercepts : ndarray, default=np.nan
        Monthly intercept values for each basin.
    monthlyvals_trends : ndarray, default=np.nan
        Monthly trend values for each basin.
    monthlyvals_numbreaks : int, default=0
        Number of breakpoints for monthly values.
    monthlyvals_datebreaks : ndarray, default=np.nan
        Dates at which the monthly value breakpoints occur.
    basin_id : ndarray, default=np.nan
        Basin number assigned to each element.
    subglacial_discharge : ndarray, default=np.nan
        Sum of subglacial discharge for each basin [m/d].
    isdischargearma : int, default=0
        Whether an ARMA model is also used for the subglacial discharge (if 0: subglacial_discharge is used, if 1: sd_ parameters are used).
    sd_ar_order : int, default=0
        Order of the subglacial discharge autoregressive model.
    sd_ma_order : int, default=0
        Order of the subglacial discharge moving-average model.
    sd_arma_timestep : int, default=0
        Time resolution of the subglacial discharge ARMA model [yr].
    sd_arlag_coefs : ndarray, default=np.nan
        Basin-specific vectors of AR lag coefficients for subglacial discharge.
    sd_malag_coefs : ndarray, default=np.nan
        Basin-specific vectors of MA lag coefficients for subglacial discharge.
    sd_monthlyfrac : ndarray, default=np.nan
        Basin-specific vectors of 12 values with fraction of the annual discharge occurring every month.
    sd_num_breaks : int, default=0
        Number of different breakpoints in the subglacial discharge piecewise-polynomial (separating sd_num_breaks+1 periods).
    sd_num_params : int, default=0
        Number of different parameters in the subglacial discharge piecewise-polynomial.
    sd_polynomialparams : ndarray, default=np.nan
        Coefficients for the subglacial discharge polynomial (const, trend, quadratic, etc.).
    sd_datebreaks : ndarray, default=np.nan
        Dates at which the breakpoints in the subglacial discharge piecewise polynomial occur (1 row per basin) [yr].

    Methods
    -------
    __init__(self, other=None)
        Initializes the RignotARMA frontalforcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the RignotARMA frontalforcings parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.frontalforcings = pyissm.build.frontalforcings.rignotarma()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.num_basins = 0
        self.num_params = 0
        self.num_breaks = 0
        self.polynomialparams = np.nan
        self.datebreaks       = np.nan
        self.ar_order = 0
        self.ma_order = 0
        self.arma_timestep = 0
        self.arlag_coefs = np.nan
        self.malag_coefs = np.nan
        self.monthlyvals_intercepts = np.nan
        self.monthlyvals_trends = np.nan
        self.monthlyvals_numbreaks = 0
        self.monthlyvals_datebreaks = np.nan
        self.basin_id = np.nan
        self.subglacial_discharge = np.nan
        self.isdischargearma = 0
        self.sd_ar_order = 0.
        self.sd_ma_order = 0.
        self.sd_arma_timestep = 0
        self.sd_arlag_coefs = np.nan
        self.sd_malag_coefs = np.nan
        self.sd_monthlyfrac = np.nan
        self.sd_num_breaks  = 0
        self.sd_num_params  = 0
        self.sd_polynomialparams = np.nan
        self.sd_datebreaks = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Frontalforcings parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'num_basins', 'number of different basins [unitless]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'basin_id', 'basin number assigned to each element [unitless]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'num_breaks', 'number of different breakpoints in the piecewise-polynomial (separating num_breaks+1 periods)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'num_params', 'number of different parameters in the piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'polynomialparams', 'coefficients for the polynomial (const,trend,quadratic,etc.),dim1 for basins,dim2 for periods,dim3 for orders, ex: polyparams=cat(num_params,intercepts,trendlinearcoefs,trendquadraticcoefs)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'datebreaks', 'dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'ar_order', 'order of the autoregressive model [unitless]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'ma_order', 'order of the moving-average model [unitless]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'arma_timestep', 'time resolution of the ARMA model [yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'arlag_coefs', 'basin-specific vectors of AR lag coefficients [unitless]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'malag_coefs', 'basin-specific vectors of MA lag coefficients [unitless]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isdischargearma','whether an ARMA model is also used for the subglacial discharge (if 0: subglacial_discharge is used, if 1: sd_ parameters are used)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'subglacial_discharge', 'sum of subglacial discharge for each basin [m/d]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sd_ar_order','order of the subglacial discharge autoregressive model [unitless]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sd_ma_order','order of the subglacial discharge moving-average model [unitless]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sd_arma_timestep','time resolution of the subglacial discharge autoregressive model [yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sd_arlag_coefs','basin-specific vectors of AR lag coefficients for subglacial discharge [unitless]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sd_malag_coefs','basin-specific vectors of MA lag coefficients for subglacial discharge [unitless]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sd_monthlyfrac','basin-specific vectors of 12 values with fraction of the annual discharge occuring every month [unitless]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sd_num_params','number of different parameters in the subglacial discharge piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sd_num_breaks','number of different breakpoints in the subglacial discharge piecewise-polynomial (separating sd_num_breaks+1 periods)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sd_datebreaks','dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sd_polynomialparams','coefficients for the sd_polynomial (const,trend,quadratic,etc.),dim1 for basins,dim2 for periods,dim3 for orders'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - frontalforcings.rignotarma Class'
        return s