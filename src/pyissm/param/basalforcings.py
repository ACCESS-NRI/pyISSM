import numpy as np
from . import param_utils
from . import class_registry

## ------------------------------------------------------
## basalforcings.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    """
    Default basal forcings parameters class for ISSM.

    This class encapsulates the default parameters for basal forcings in the ISSM (Ice Sheet System Model) framework.
    It defines the melting rates for grounded and floating ice, perturbation melting rate, and geothermal heat flux.
    
    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    groundedice_melting_rate : ndarray, default=np.nan
        Basal melting rate for grounded ice (positive if melting) [m/yr].
    floatingice_melting_rate : ndarray, default=np.nan
        Basal melting rate for floating ice (positive if melting) [m/yr].
    perturbation_melting_rate : ndarray, default=np.nan
        Optional perturbation in basal melting rate under floating ice (positive if melting) [m/yr].
    geothermalflux : float, default=np.nan
        Geothermal heat flux [W/m^2].

    Methods
    -------
    __init__(self, other=None)
        Initializes the basal forcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the basal forcings parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.basalforcings = pyissm.build.basalforcings.default()
    md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices,))
    md.basalforcings.floatingice_melting_rate = np.ones((md.mesh.numberofvertices,)) * 2
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.groundedice_melting_rate = np.nan
        self.floatingice_melting_rate = np.nan
        self.perturbation_melting_rate = np.nan
        self.geothermalflux = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   basal forcings parameters:\n'
        
        s += '{}\n'.format(param_utils.fielddisplay(self, 'groundedice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'floatingice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'perturbation_melting_rate', '(optional) perturbation in basal melting rate under floating ice [m/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'geothermalflux', 'geothermal heat flux [W/m^2]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - basalforcings.default Class'
        return s

## ------------------------------------------------------
## basalforcings.pico
## ------------------------------------------------------
@class_registry.register_class
class pico(class_registry.manage_state):
    """
    Potsdam Ice-shelf Cavity mOdel (PICO) basal forcings parameters class for ISSM.

    This class encapsulates the parameters for the PICO basal melt parameterization in the ISSM (Ice Sheet System Model) framework.
    It defines the structure of the ice shelf cavities, including the number of basins, basin IDs, and various parameters related to ocean temperature, salinity, and melting rates.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    num_basins : int, default=0
        Number of basins the model domain is partitioned into [unitless].
    basin_id : ndarray, default=np.nan
        Basin number assigned to each node [unitless].
    maxboxcount : int, default=0
        Maximum number of boxes initialized under all ice shelves.
    overturning_coeff : float, default=np.nan
        Overturning strength [m^3/s].
    gamma_T : float, default=0.
        Turbulent temperature exchange velocity [m/s].
    farocean_temperature : ndarray, default=np.nan
        Depth averaged ocean temperature in front of the ice shelf for each basin [K].
    farocean_salinity : ndarray, default=np.nan
        Depth averaged ocean salinity in front of the ice shelf for each basin [psu].
    isplume : int, default=0
        Boolean (0 or 1) to use buoyant plume melt rate parameterization from Lazeroms et al., 2018 (default false).
    geothermalflux : float, default=np.nan
        Geothermal heat flux [W/m^2].
    groundedice_melting_rate : ndarray, default=np.nan
        Basal melting rate for grounded ice (positive if melting) [m/yr].

    Methods
    -------
    __init__(self, other=None)
        Initializes the PICO basal forcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the PICO basal forcings parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.basalforcings = pyissm.build.basalforcings.pico()
    md.basalforcings.num_basins = 3
    md.basalforcings.basin_id = np.array([1, 2, 3])
    md.basalforcings.farocean_temperature = np.array([273.15, 273.2, 273.1])
    """

    # Initialise with default parameters
    def __init__(self, other=None):
        self.num_basins = 0
        self.basin_id = np.nan
        self.maxboxcount = 0
        self.overturning_coeff = np.nan
        self.gamma_T = 0.
        self.farocean_temperature = np.nan
        self.farocean_salinity = np.nan
        self.isplume = 0
        self.geothermalflux = np.nan
        self.groundedice_melting_rate = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   PICO basal melt rate parameterization:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self,'num_basins','number of basins the model domain is partitioned into [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'basin_id','basin number assigned to each node [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'maxboxcount','maximum number of boxes initialized under all ice shelves'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'overturning_coeff','overturning strength [m^3/s]'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'gamma_T','turbulent temperature exchange velocity [m/s]'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'farocean_temperature','depth averaged ocean temperature in front of the ice shelf for basin i [K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'farocean_salinity','depth averaged ocean salinity in front of the ice shelf for basin i [psu]'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'isplume','boolean to use buoyant plume melt rate parameterization from Lazeroms et al., 2018 (default false)'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'geothermalflux','geothermal heat flux [W/m^2]'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'groundedice_melting_rate','basal melting rate (positive if melting) [m/yr]'))

        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - basalforcings.pico Class'
        return s

## ------------------------------------------------------
## basalforcings.linear
## ------------------------------------------------------
@class_registry.register_class
class linear(class_registry.manage_state):
    """
    Linear basal forcings parameters class for ISSM.

    This class encapsulates the parameters for linear basal forcings in the ISSM (Ice Sheet System Model) framework.
    It defines the melting rates for deep and upper water, grounded ice, and geothermal flux, allowing for a depth-dependent representation of basal melting processes.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    deepwater_melting_rate : float, default=0.
        Basal melting rate applied for floating ice with base < deepwater_elevation [m/yr].
    deepwater_elevation : float, default=0.
        Elevation threshold for deepwater melting rate [m].
    upperwater_melting_rate : float, default=0.
        Basal melting rate applied for floating ice with base >= upperwater_elevation [m/yr].
    upperwater_elevation : float, default=0.
        Elevation threshold for upperwater melting rate [m].
    groundedice_melting_rate : ndarray, default=np.nan
        Basal melting rate for grounded ice (positive if melting) [m/yr].
    perturbation_melting_rate : ndarray, default=np.nan
        Perturbation applied to computed melting rate (positive if melting) [m/yr].
    geothermalflux : float, default=np.nan
        Geothermal heat flux [W/m^2].

    Methods
    -------
    __init__(self, other=None)
        Initializes the linear basal forcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the linear basal forcings parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.basalforcings = pyissm.build.basalforcings.linear()
    md.basalforcings.deepwater_melting_rate = 1.5
    md.basalforcings.deepwater_elevation = -500
    md.basalforcings.upperwater_melting_rate = 0.5
    md.basalforcings.upperwater_elevation = -200
    """

    # Initialise with default parameters
    def __init__(self, other=None):
        self.deepwater_melting_rate = 0.
        self.deepwater_elevation = 0.
        self.upperwater_melting_rate = 0.
        self.upperwater_elevation = 0.
        self.groundedice_melting_rate = np.nan
        self.perturbation_melting_rate = np.nan
        self.geothermalflux = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   linear basal forcings parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, "deepwater_melting_rate", "basal melting rate (positive if melting applied for floating ice whith base < deepwater_elevation) [m/yr]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "deepwater_elevation", "elevation of ocean deepwater [m]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "upperwater_melting_rate", "upper melting rate (positive if melting applied for floating ice whith base >= upperwater_elevation) [m/yr]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "upperwater_elevation", "elevation of ocean upper water [m]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "groundedice_melting_rate", "basal melting rate (positive if melting) [m/yr]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "perturbation_melting_rate", "perturbation applied to computed melting rate (positive if melting) [m/yr]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "geothermalflux", "geothermal heat flux [W/m^2]"))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - basalforcings.linear Class'
        return s

## ------------------------------------------------------
## basalforcings.lineararma
## ------------------------------------------------------
@class_registry.register_class
class lineararma(class_registry.manage_state):
    """
    Linear ARMA basal forcings parameters class for ISSM.

    This class encapsulates the parameters for linear ARMA basal forcings in the ISSM (Ice Sheet System Model) framework.
    It defines the structure for piecewise polynomial parameters, autoregressive and moving-average coefficients, and various melting rates and elevations for different basins.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    num_basins : int, default=0
        Number of different basins [unitless].
    num_params : int, default=0
        Number of different parameters in the piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.).
    num_breaks : int, default=0
        Number of different breakpoints in the piecewise-polynomial (separating num_breaks+1 periods).
    polynomialparams : ndarray, default=np.nan
        Coefficients for the polynomial (const, trend, quadratic, etc.), dim1 for basins, dim2 for periods, dim3 for orders.
    datebreaks : ndarray, default=np.nan
        Dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr].
    ar_order : float, default=0.
        Order of the autoregressive model [unitless].
    ma_order : float, default=0.
        Order of the moving-average model [unitless].
    arma_timestep : int, default=0
        Time resolution of the ARMA model [yr].
    arlag_coefs : ndarray, default=np.nan
        Basin-specific vectors of AR lag coefficients [unitless].
    malag_coefs : ndarray, default=np.nan
        Basin-specific vectors of MA lag coefficients [unitless].
    basin_id : ndarray, default=np.nan
        Basin number assigned to each element [unitless].
    groundedice_melting_rate : ndarray, default=np.nan
        Node-specific basal melting rate for grounded ice (positive if melting) [m/yr].
    deepwater_elevation : ndarray, default=np.nan
        Basin-specific elevation of ocean deepwater [m].
    upperwater_melting_rate : ndarray, default=np.nan
        Basin-specific basal melting rate (positive if melting applied for floating ice with base >= upperwater_elevation) [m/yr].
    upperwater_elevation : ndarray, default=np.nan
        Basin-specific elevation of ocean upperwater [m].
    geothermalflux : ndarray, default=np.nan
        Node-specific geothermal heat flux [W/m^2].

    Methods
    -------
    __init__(self, other=None)
        Initializes the linear ARMA basal forcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the linear ARMA basal forcings parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.basalforcings = pyissm.build.basalforcings.lineararma()
    """

    # Initialise with default parameters
    def __init__(self, other=None):
        self.num_basins = 0
        self.num_params = 0
        self.num_breaks = 0
        self.polynomialparams = np.nan
        self.datebreaks       = np.nan
        self.ar_order = 0.
        self.ma_order = 0.
        self.arma_timestep = 0
        self.arlag_coefs = np.nan
        self.malag_coefs = np.nan
        self.basin_id = np.nan
        self.groundedice_melting_rate = np.nan
        self.deepwater_elevation = np.nan
        self.upperwater_melting_rate = np.nan
        self.upperwater_elevation = np.nan
        self.geothermalflux = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   basal forcings parameters:\n'

        s += '   autoregressive model is applied for deepwater_melting_rate\n'
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
        s += '{}\n'.format(param_utils.fielddisplay(self, 'deepwater_elevation', 'basin-specific elevation of ocean deepwater [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'upperwater_melting_rate', 'basin-specic basal melting rate (positive if melting applied for floating ice whith base >= upperwater_elevation) [m/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'upperwater_elevation', 'basin-specific elevation of ocean upperwater [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'groundedice_melting_rate','node-specific basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'geothermalflux','node-specific geothermal heat flux [W/m^2]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - basalforcings.lineararmaClass'
        return s

## ------------------------------------------------------
## basalforcings.mismip
## ------------------------------------------------------
@class_registry.register_class
class mismip(class_registry.manage_state):
    """
    MISMIP basal forcings parameters class for ISSM.

    This class encapsulates the parameters for the MISMIP basal melt parameterization in the ISSM (Ice Sheet System Model) framework.
    It defines the basal melting rate for grounded ice, a melt rate factor, a threshold thickness for saturation of basal melting,
    an upper depth above which the melt rate is zero, and geothermal heat flux.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    groundedice_melting_rate : ndarray, default=np.nan
        Basal melting rate for grounded ice (positive if melting) [m/yr].
    meltrate_factor : float, default=0.2
        Melt-rate factor [1/yr] (sign is opposite to MISMIP+ benchmark to remain consistent with ISSM convention of positive values for melting).
    threshold_thickness : float, default=75.
        Threshold thickness for saturation of basal melting [m].
    upperdepth_melt : float, default=-100.
        Depth above which melt rate is zero [m].
    geothermalflux : float, default=np.nan
        Geothermal heat flux [W/m^2].

    Methods
    -------
    __init__(self, other=None)
        Initializes the MISMIP basal forcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the MISMIP basal forcings parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.basalforcings = pyissm.build.basalforcings.mismip()
    md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices,))
    md.basalforcings.meltrate_factor = 0.2
    md.basalforcings.threshold_thickness = 75.
    md.basalforcings.upperdepth_melt = -100.
    """

    # Initialise with default parameters
    def __init__(self, other=None):
        self.groundedice_melting_rate = np.nan
        self.meltrate_factor = 0.2
        self.threshold_thickness = 75.
        self.upperdepth_melt = -100.
        self.geothermalflux = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   MISMIP + basal melt parameterization\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, "groundedice_melting_rate", "basal melting rate (positive if melting) [m / yr]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "meltrate_factor", "Melt - rate rate factor [1 / yr] (sign is opposite to MISMIP + benchmark to remain consistent with ISSM convention of positive values for melting)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "threshold_thickness", "Threshold thickness for saturation of basal melting [m]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "upperdepth_melt", "Depth above which melt rate is zero [m]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "geothermalflux", "Geothermal heat flux [W / m^2]"))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - basalforcings.mismip Class'
        return s

## ------------------------------------------------------
## basalforcings.plume
## ------------------------------------------------------
@class_registry.register_class
class plume(class_registry.manage_state):
    """
    Plume basal forcings parameters class for ISSM.

    This class encapsulates the parameters for plume basal forcings in the ISSM (Ice Sheet System Model) framework.
    It defines the structure of the mantle plume, including its radius, depth, and position, as well as parameters related to geothermal heat flux and melting rates.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    floatingice_melting_rate : ndarray, default=np.nan
        Basal melting rate for floating ice (positive if melting) [m/yr].
    groundedice_melting_rate : ndarray, default=np.nan
        Basal melting rate for grounded ice (positive if melting) [m/yr].
    mantleconductivity : float, default=2.2
        Mantle heat conductivity [W/m^3].
    nusselt : float, default=300
        Nusselt number, ratio of mantle to plume [1].
    dtbg : float, default=0.011
        Background temperature gradient [degree/m].
    plumeradius : float, default=100000
        Radius of the mantle plume [m].
    topplumedepth : float, default=10000
        Depth of the mantle plume top below the crust [m].
    bottomplumedepth : float, default=1050000
        Depth of the mantle plume base below the crust [m].
    plumex : float, default=np.nan
        x coordinate of the center of the plume [m].
    plumey : float, default=np.nan
        y coordinate of the center of the plume [m].
    crustthickness : float, default=30000
        Thickness of the crust [m].
    uppercrustthickness : float, default=14000
        Thickness of the upper crust [m].
    uppercrustheat : float, default=1.7e-6
        Volumic heat of the upper crust [W/m^3].
    lowercrustheat : float, default=0.4e-6
        Volumic heat of the lower crust [W/m^3].

    Methods
    -------
    __init__(self, other=None)
        Initializes the plume basal forcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the plume basal forcings parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.basalforcings = pyissm.build.basalforcings.plume()
    md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices,))
    md.basalforcings.floatingice_melting_rate = np.ones((md.mesh.numberofvertices,)) * 2
    """

    # Initialise with default parameters
    def __init__(self, other=None):
        self.floatingice_melting_rate = np.nan
        self.groundedice_melting_rate = np.nan
        self.mantleconductivity = 2.2
        self.nusselt = 300
        self.dtbg = 11 / 1000.
        self.plumeradius = 100000
        self.topplumedepth = 10000
        self.bottomplumedepth = 1050000
        self.plumex = np.nan
        self.plumey = np.nan
        self.crustthickness = 30000
        self.uppercrustthickness = 14000
        self.uppercrustheat = 1.7 * pow(10, -6)
        self.lowercrustheat = 0.4 * pow(10, -6)

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   mantle plume basal melt parameterization:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'groundedice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'floatingice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'mantleconductivity', 'mantle heat conductivity [W/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'nusselt', 'nusselt number, ratio of mantle to plume [1]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'dtbg', 'background temperature gradient [degree/m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'plumeradius', 'radius of the mantle plume [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'topplumedepth', 'depth of the mantle plume top below the crust [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'bottomplumedepth', 'depth of the mantle plume base below the crust [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'plumex', 'x coordinate of the center of the plume [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'plumey', 'y coordinate of the center of the plume [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'crustthickness', 'thickness of the crust [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'uppercrustthickness', 'thickness of the upper crust [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'uppercrustheat', 'volumic heat of the upper crust [w/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'lowercrustheat', 'volumic heat of the lowercrust [w/m^3]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - basalforcings.plume Class'
        return s

## ------------------------------------------------------
## basalforcings.spatiallinear
## ------------------------------------------------------
@class_registry.register_class
class spatiallinear(class_registry.manage_state):
    """
    Spatial linear basal forcings parameters class for ISSM.

    This class encapsulates the parameters for spatial linear basal forcings in the ISSM (Ice Sheet System Model) framework.
    It defines the melting rates for grounded ice, deepwater, and upperwater, as well as geothermal heat flux and perturbation melting rate.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    groundedice_melting_rate : ndarray, default=np.nan
        Basal melting rate for grounded ice (positive if melting) [m/yr].
    deepwater_melting_rate : ndarray, default=np.nan
        Basal melting rate applied for floating ice with base < deepwater_elevation [m/yr].
    deepwater_elevation : ndarray, default=np.nan
        Elevation threshold for deepwater melting rate [m].
    upperwater_melting_rate : ndarray, default=np.nan
        Basal melting rate applied for floating ice with base >= upperwater_elevation [m/yr].
    upperwater_elevation : ndarray, default=np.nan
        Elevation threshold for upperwater melting rate [m].
    geothermalflux : ndarray, default=np.nan
        Geothermal heat flux [W/m^2].
    perturbation_melting_rate : ndarray, default=np.nan
        Basal melting rate perturbation added to computed melting rate (positive if melting) [m/yr].

    Methods
    -------
    __init__(self, other=None)
        Initializes the spatial linear basal forcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the spatial linear basal forcings parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.basalforcings = pyissm.build.basalforcings.spatiallinear()
    """

    # Initialise with default parameters
    def __init__(self, other=None):
        self.groundedice_melting_rate = np.nan
        self.deepwater_melting_rate = np.nan
        self.deepwater_elevation = np.nan
        self.upperwater_melting_rate = np.nan
        self.upperwater_elevation = np.nan
        self.geothermalflux = np.nan
        self.perturbation_melting_rate = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   spatial linear basal forcings parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'groundedice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'deepwater_melting_rate', 'basal melting rate (positive if melting applied for floating ice whith base < deepwater_elevation) [m/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'deepwater_elevation', 'elevation of ocean deepwater [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'upperwater_melting_rate', 'basal melting rate (positive if melting applied for floating ice whith base >= upperwater_elevation) [m/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'upperwater_elevation', 'elevation of ocean upperwater [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'perturbation_melting_rate', 'basal melting rate perturbation added to computed melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'geothermalflux', 'geothermal heat flux [W/m^2]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - basalforcings.spatiallinear Class'
        return s
