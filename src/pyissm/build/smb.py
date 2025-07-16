import numpy as np
from . import build_utils
from . import class_registry

## ------------------------------------------------------
## smb.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    '''
    smb.default Class definition
    '''

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

        s += '{}\n'.format(build_utils.fielddisplay(self, 'mass_balance', 'surface mass balance [m/yr ice eq]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
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
    '''
    smb.arma Class definition
    '''

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
        s += '{}\n'.format(build_utils.fielddisplay(self, 'lapserates', 'basin-specific SMB lapse rates applied in each elevation bin, 1 row per basin, 1 column per bin, dimension 3 can be of size 12 to prescribe monthly varying values [m ice eq yr^-1 m^-1] (default: no lapse rate)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'elevationbins', 'basin-specific separations between elevation bins, 1 row per basin, 1 column per limit between bins, dimension 3 can be of size 12 to prescribe monthly varying values [m] (default: no basin separation)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'refelevation', 'basin-specific reference elevations at which SMB is calculated, and from which SMB is downscaled using lapserates (default: basin mean elevation) [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
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
    '''
    smb.components Class definition
    '''

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

        s += '{}\n'.format(build_utils.fielddisplay(self, 'accumulation', 'accumulated snow [m/yr ice eq]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'runoff', 'amount of ice melt lost from the ice column [m/yr ice eq]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'evaporation', 'mount of ice lost to evaporative processes [m/yr ice eq]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
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
    '''
    smb.d18opdd Class definition
    '''

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

        s += '{}\n'.format(build_utils.fielddisplay(self, 'isd18opd', 'is delta18o parametrisation from present day temperature and precipitation activated (0 or 1, default is 0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'issetpddfac', 'is user passing in defined pdd factors (0 or 1, default is 0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'desfac', 'desertification elevation factor (between 0 and 1, default is 0.5) [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 's0p', 'should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 's0t', 'should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rlaps', 'present day lapse rate [degree/km]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'temperatures_presentday', 'monthly present day surface temperatures [K], required if delta18o/mungsm is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'precipitations_presentday', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'istemperaturescaled', 'if delta18o parametrisation from present day temperature and precipitation is activated, is temperature scaled to delta18o value? (0 or 1, default is 1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isprecipscaled', 'if delta18o parametrisation from present day temperature and precipitation is activated, is precipitation scaled to delta18o value? (0 or 1, default is 1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'temperatures_reconstructed', 'monthly historical surface temperatures [K], required if delta18o/mungsm/d18opd is activated and istemperaturescaled is not activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'precipitations_reconstructed', 'monthly historical precipitation [m/yr water eq], required if delta18o/mungsm/d18opd is activated and isprecipscaled is not activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'delta18o', 'delta18o [per mil], required if pdd is activated and delta18o activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'dpermil', 'degree per mil, required if d18opd is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'f', 'precip/temperature scaling factor, required if d18opd is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'pddfac_snow', 'Pdd factor for snow for all the domain [mm ice equiv/day/degree C]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'pddfac_ice', 'Pdd factor for ice for all the domain [mm ice equiv/day/degree C]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
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
    '''
    smb.gradients Class definition
    '''

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

        s += '{}\n'.format(build_utils.fielddisplay(self, 'issmbgradients', 'is smb gradients method activated (0 or 1, default is 0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'href', 'reference elevation from which deviation is used to calculate SMB adjustment in smb gradients method'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'smbref', 'reference smb from which deviation is calculated in smb gradients method [m/yr ice equiv]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'b_pos', 'slope of hs - smb regression line for accumulation regime required if smb gradients is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'b_neg', 'slope of hs - smb regression line for ablation regime required if smb gradients is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
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
    '''
    smb.gradientscomponents Class definition
    '''

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

        s += '{}\n'.format(build_utils.fielddisplay(self, 'issmbgradients', 'is smb gradients method activated (0 or 1, default is 0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'accuref', ' reference value of the accumulation'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'accualti', ' Altitude at which the accumulation is equal to the reference value'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'accugrad', ' Gradient of the variation of the accumulation (0 for uniform accumulation)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'runoffref', ' reference value of the runoff m w.e. y-1 (temperature times ddf)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'runoffalti', ' Altitude at which the runoff is equal to the reference value'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'runoffgrad', ' Gradient of the variation of the runoff (0 for uniform runoff) m w.e. m-1 y-1 (lapse rate times ddf)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
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
    '''
    smb.gradientsela Class definition
    '''

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
        s += '{}\n'.format(build_utils.fielddisplay(self, 'ela', ' equilibrium line altitude from which deviation is used to calculate smb using the smb gradients ela method [m a.s.l.]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'b_pos', ' vertical smb gradient (dB/dz) above ela'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'b_neg', ' vertical smb gradient (dB/dz) below ela'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'b_max', ' upper cap on smb rate, default: 9999 (no cap) [m ice eq./yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'b_min', ' lower cap on smb rate, default: -9999 (no cap) [m ice eq./yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
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
    '''
    smb.henning Class definition
    '''

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

        s += '{}\n'.format(build_utils.fielddisplay(self, 'mass_balance', 'surface mass balance [m/yr ice eq]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
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
    '''
    smb.meltcomponents Class definition
    '''

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

        s += '{}\n'.format(build_utils.fielddisplay(self, 'accumulation', 'accumulated snow [m/yr ice eq]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'evaporation', 'mount of ice lost to evaporative processes [m/yr ice eq]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'melt', 'amount of ice melt in the ice column [m/yr ice eq]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'refreeze', 'amount of ice melt refrozen in the ice column [m/yr ice eq]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
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
    '''
    smb.pdd Class definition
    '''

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

        s += '{}\n'.format(build_utils.fielddisplay(self, 'isdelta18o', 'is temperature and precipitation delta18o parametrisation activated (0 or 1, default is 0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'ismungsm', 'is temperature and precipitation mungsm parametrisation activated (0 or 1, default is 0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'issetpddfac', 'is user passing in defined pdd factors (0 or 1, default is 0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'desfac', 'desertification elevation factor (between 0 and 1, default is 0.5) [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 's0p', 'should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 's0t', 'should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rlaps', 'present day lapse rate [degree/km]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rlapslgm', 'LGM lapse rate [degree/km]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'monthlytemperatures', 'monthly surface temperatures [K], required if pdd is activated and delta18o not activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'precipitation', 'monthly surface precipitation [m/yr water eq], required if pdd is activated and delta18o or mungsm not activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'delta18o', 'delta18o [per mil], required if pdd is activated and delta18o activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'delta18o_surface', 'surface elevation of the delta18o site, required if pdd is activated and delta18o activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'temperatures_presentday', 'monthly present day surface temperatures [K], required if delta18o/mungsm is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'temperatures_lgm', 'monthly LGM surface temperatures [K], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'precipitations_presentday', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'precipitations_lgm', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'Tdiff', 'time interpolation parameter for temperature, 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sealev', 'sea level [m], 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'temperatures_presentday', 'monthly present day surface temperatures [K], required if delta18o/mungsm is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'temperatures_lgm', 'monthly LGM surface temperatures [K], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'precipitations_presentday', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'precipitations_lgm', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'Pfac', 'time interpolation parameter for precipitation, 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'Tdiff', 'time interpolation parameter for temperature, 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sealev', 'sea level [m], 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
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
    '''
    smb.pddSicopolis Class definition
    '''

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
        s += '{}\n'.format(build_utils.fielddisplay(self, 'monthlytemperatures', 'monthly surface temperatures [K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'precipitation', 'monthly surface precipitation [m/yr water eq]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'temperature_anomaly', 'anomaly to monthly reference temperature (additive [K])'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'precipitation_anomaly', 'anomaly to monthly precipitation (multiplicative, e.g. q = q0*exp(0.070458*DeltaT) after Huybrechts (2002)) [no unit])'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'smb_corr', 'correction of smb after PDD call [m/a]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 's0p', 'should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 's0t', 'should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rlaps', 'present day lapse rate (default is 7.4 degree/km)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'desfac', 'desertification elevation factor (default is -log(2.0)/1000)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isfirnwarming', 'is firnwarming (Reeh 1991) activated (0 or 1, default is 1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested (TemperaturePDD, SmbAccumulation, SmbMelt)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.pddSicopolis Class'
        return s