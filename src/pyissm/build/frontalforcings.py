import numpy as np
from . import build_utils
from . import class_registry

## ------------------------------------------------------
## frontalforcings.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    '''
    frontalforcings.default Class definition
    '''

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
    '''
    frontalforcings.rignot Class definition
    '''

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
    '''
    frontalforcings.rignotarma Class definition
    '''

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