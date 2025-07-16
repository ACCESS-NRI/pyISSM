import numpy as np
from . import build_utils
from . import class_registry

## ------------------------------------------------------
## basalforcings.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    '''
    basalforcings.default Class definition
    '''

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
        
        s += '{}\n'.format(build_utils.fielddisplay(self, 'groundedice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'floatingice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'perturbation_melting_rate', '(optional) perturbation in basal melting rate under floating ice [m/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'geothermalflux', 'geothermal heat flux [W/m^2]'))
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
    '''
    basalforcings.pico Class definition
    '''

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

        s += '{}\n'.format(build_utils.fielddisplay(self,'num_basins','number of basins the model domain is partitioned into [unitless]'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'basin_id','basin number assigned to each node [unitless]'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'maxboxcount','maximum number of boxes initialized under all ice shelves'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'overturning_coeff','overturning strength [m^3/s]'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'gamma_T','turbulent temperature exchange velocity [m/s]'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'farocean_temperature','depth averaged ocean temperature in front of the ice shelf for basin i [K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'farocean_salinity','depth averaged ocean salinity in front of the ice shelf for basin i [psu]'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'isplume','boolean to use buoyant plume melt rate parameterization from Lazeroms et al., 2018 (default false)'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'geothermalflux','geothermal heat flux [W/m^2]'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'groundedice_melting_rate','basal melting rate (positive if melting) [m/yr]'))

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
    '''
    basalforcings.linear Class definition
    '''

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

        s += '{}\n'.format(build_utils.fielddisplay(self, "deepwater_melting_rate", "basal melting rate (positive if melting applied for floating ice whith base < deepwater_elevation) [m/yr]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "deepwater_elevation", "elevation of ocean deepwater [m]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "upperwater_melting_rate", "upper melting rate (positive if melting applied for floating ice whith base >= upperwater_elevation) [m/yr]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "upperwater_elevation", "elevation of ocean upper water [m]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "groundedice_melting_rate", "basal melting rate (positive if melting) [m/yr]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "perturbation_melting_rate", "perturbation applied to computed melting rate (positive if melting) [m/yr]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "geothermalflux", "geothermal heat flux [W/m^2]"))
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
    '''
    basalforcings.lineararma Class definition
    '''

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
        s += '{}\n'.format(build_utils.fielddisplay(self, 'deepwater_elevation', 'basin-specific elevation of ocean deepwater [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'upperwater_melting_rate', 'basin-specic basal melting rate (positive if melting applied for floating ice whith base >= upperwater_elevation) [m/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'upperwater_elevation', 'basin-specific elevation of ocean upperwater [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'groundedice_melting_rate','node-specific basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'geothermalflux','node-specific geothermal heat flux [W/m^2]'))
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
    '''
    basalforcings.mismip Class definition
    '''

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

        s += '{}\n'.format(build_utils.fielddisplay(self, "groundedice_melting_rate", "basal melting rate (positive if melting) [m / yr]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "meltrate_factor", "Melt - rate rate factor [1 / yr] (sign is opposite to MISMIP + benchmark to remain consistent with ISSM convention of positive values for melting)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "threshold_thickness", "Threshold thickness for saturation of basal melting [m]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "upperdepth_melt", "Depth above which melt rate is zero [m]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "geothermalflux", "Geothermal heat flux [W / m^2]"))
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
    '''
    basalforcings.plume Class definition
    '''

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

        s += '{}\n'.format(build_utils.fielddisplay(self, 'groundedice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'floatingice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'mantleconductivity', 'mantle heat conductivity [W/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'nusselt', 'nusselt number, ratio of mantle to plume [1]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'dtbg', 'background temperature gradient [degree/m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'plumeradius', 'radius of the mantle plume [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'topplumedepth', 'depth of the mantle plume top below the crust [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'bottomplumedepth', 'depth of the mantle plume base below the crust [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'plumex', 'x coordinate of the center of the plume [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'plumey', 'y coordinate of the center of the plume [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'crustthickness', 'thickness of the crust [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'uppercrustthickness', 'thickness of the upper crust [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'uppercrustheat', 'volumic heat of the upper crust [w/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'lowercrustheat', 'volumic heat of the lowercrust [w/m^3]'))
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
    '''
    basalforcings.spatiallinear Class definition
    '''

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

        s += '{}\n'.format(build_utils.fielddisplay(self, 'groundedice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'deepwater_melting_rate', 'basal melting rate (positive if melting applied for floating ice whith base < deepwater_elevation) [m/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'deepwater_elevation', 'elevation of ocean deepwater [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'upperwater_melting_rate', 'basal melting rate (positive if melting applied for floating ice whith base >= upperwater_elevation) [m/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'upperwater_elevation', 'elevation of ocean upperwater [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'perturbation_melting_rate', 'basal melting rate perturbation added to computed melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'geothermalflux', 'geothermal heat flux [W/m^2]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - basalforcings.spatiallinear Class'
        return s
