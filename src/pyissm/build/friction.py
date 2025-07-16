import numpy as np
from . import build_utils
from . import class_registry

## ------------------------------------------------------
## friction.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    '''
    friction.default Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.coefficient = np.nan
        self.p = np.nan
        self.q = np.nan
        self.coupling = 0
        self.linearize = 0
        self.effective_pressure = np.nan
        self.effective_pressure_limit = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Basal shear stress parameters: Sigma_b = coefficient^2 * Neff ^r * |u_b|^(s - 1) * u_b,\n'
        s += '(effective stress Neff = rho_ice * g * thickness + rho_water * g * base, r = q / p and s = 1 / p)\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'coefficient', 'friction coefficient [SI]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'p', 'p exponent'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'q', 'q exponent'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'coupling', 'Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: used coupled model (not implemented yet)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'linearize', '0: not linearized, 1: interpolated linearly, 2: constant per element (default is 0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.default Class'
        return s

## ------------------------------------------------------
## friction.coulomb
## ------------------------------------------------------
@class_registry.register_class
class coulomb(class_registry.manage_state):
    '''
    friction.coulomb Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.coefficient = np.nan
        self.coefficientcoulomb = np.nan
        self.p = np.nan
        self.q = np.nan
        self.coupling = 0
        self.effective_pressure = np.nan
        self.effective_pressure_limit = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Basal shear stress parameters: Sigma_b = min(coefficient^2 * Neff ^r * |u_b|^(s - 1) * u_b,\n'
        s += 'coefficientcoulomb^2 * Neff), (effective stress Neff = rho_ice * g * thickness + rho_water * g * bed, r = q / p and s = 1 / p).\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'coefficient', 'power law (Weertman) friction coefficient [SI]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'coefficientcoulomb', 'Coulomb friction coefficient [SI]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'p', 'p exponent'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'q', 'q exponent'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'coupling', 'Coupling flag: 0 for default, 1 for forcing(provide md.friction.effective_pressure)  and 2 for coupled(not implemented yet)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.coulomb Class'
        return s

## ------------------------------------------------------
## friction.coulomb2
## ------------------------------------------------------
@class_registry.register_class
class coulomb2(class_registry.manage_state):
    '''
    friction.coulomb2 Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.coefficient = np.nan
        self.coefficientcoulomb = np.nan
        self.p = np.nan
        self.q = np.nan
        self.coupling = 0
        self.effective_pressure = np.nan
        self.effective_pressure_limit = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Basal shear stress parameters: Sigma_b = min(coefficient^2 * Neff ^r * |u_b|^(s - 1) * u_b,\n'
        s += 'coefficientcoulomb^2 * Neff), (effective stress Neff = rho_ice * g * thickness + rho_water * g * bed, r = q / p and s = 1 / p).\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'coefficient', 'power law (Weertman) friction coefficient [SI]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'coefficientcoulomb', 'Coulomb friction coefficient [SI]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'p', 'p exponent'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'q', 'q exponent'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'coupling', 'Coupling flag: 0 for default, 1 for forcing(provide md.friction.effective_pressure)  and 2 for coupled(not implemented yet)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.coulomb2 Class'
        return s

## ------------------------------------------------------
## friction.hydro
## ------------------------------------------------------
@class_registry.register_class
class hydro(class_registry.manage_state):
    '''
    friction.hydro Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.coupling = 0
        self.q = np.nan
        self.C = np.nan
        self.As = np.nan
        self.effective_pressure = np.nan
        self.effective_pressure_limit = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Effective Pressure based friction law described in Gagliardini 2007\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'coupling', 'Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: used coupled model (not implemented yet)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'q', 'friction law exponent q >= 1'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'C', 'friction law max value (Iken bound)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'As', 'Sliding Parameter without cavitation [m Pa^ - n s^ - 1]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.hydro Class'
        return s

## ------------------------------------------------------
## friction.josh
## ------------------------------------------------------
@class_registry.register_class
class josh(class_registry.manage_state):
    '''
    friction.josh Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.coefficient = np.nan
        self.pressure_adjusted_temperature = np.nan
        self.gamma = 1.
        self.effective_pressure_limit = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Basal shear stress parameters: Sigma_b = coefficient^2 * Neff ^r * |u_b|^(s - 1) * u_b,\n'
        s += '(effective stress Neff = rho_ice * g * thickness + rho_water * g * base, r = q / p and s = 1 / p)\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'coefficient', 'friction coefficient [SI]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'pressure_adjusted_temperature', 'friction pressure_adjusted_temperature (T - Tpmp) [K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'gamma', '(T - Tpmp)/gamma [K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.josh Class'
        return s

## ------------------------------------------------------
## friction.pism
## ------------------------------------------------------
@class_registry.register_class
class pism(class_registry.manage_state):
    '''
    friction.pism Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.pseudoplasticity_exponent = 0.6
        self.threshold_speed = 100.
        self.delta = 0.02
        self.void_ratio = 0.69
        self.till_friction_angle = np.nan
        self.sediment_compressibility_coefficient = np.nan
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Basal shear stress parameters for the PISM friction law (See Aschwanden et al. 2016 for more details)\n'

        s += "{}\n".format(build_utils.fielddisplay(self, 'pseudoplasticity_exponent', 'pseudoplasticity exponent [dimensionless]'))
        s += "{}\n".format(build_utils.fielddisplay(self, 'threshold_speed', 'threshold speed [m / yr]'))
        s += "{}\n".format(build_utils.fielddisplay(self, 'delta', 'lower limit of the effective pressure, expressed as a fraction of overburden pressure [dimensionless]'))
        s += "{}\n".format(build_utils.fielddisplay(self, 'void_ratio', 'void ratio at a reference effective pressure [dimensionless]'))
        s += "{}\n".format(build_utils.fielddisplay(self, 'till_friction_angle', 'till friction angle [deg], recommended default: 30 deg'))
        s += "{}\n".format(build_utils.fielddisplay(self, 'sediment_compressibility_coefficient', 'coefficient of compressibility of the sediment [dimensionless], recommended default: 0.12'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.pism Class'
        return s

## ------------------------------------------------------
## friction.regcoulomb
## ------------------------------------------------------
@class_registry.register_class
class regcoulomb(class_registry.manage_state):
    '''
    friction.regcoulomb Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.C = np.nan
        self.u0 = 1000
        self.m = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        # See Joughin et al. 2019 (equivalent form by Matt Trevers, poster at AGU 2022) https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019GL082526
        s = 'Regularized Coulomb friction law (Joughin et al., 2019) parameters:\n'

        s += '   Regularized Coulomb friction law reads:\n'
        s += '                       C^2 |u|^(1/m)         \n'
        s += '      tau_b = -  ____________________________\n'
        s += '                     (|u|/u0 + 1)^(1/m)      \n'
        s += '\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'C', 'friction coefficient [SI]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'm', 'm exponent (set to m = 3 in original paper)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'u0', 'velocity controlling plastic limit'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.regcoulomb Class'
        return s

## ------------------------------------------------------
## friction.regcoulomb2
## ------------------------------------------------------
@class_registry.register_class
class regcoulomb2(class_registry.manage_state):
    '''
    friction.regcoulomb2 Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.C = np.nan
        self.K = np.nan
        self.m = np.nan
        self.effective_pressure_limit = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        # See Zoet and Iverson 2020 or Choi et al., 2022
        s = 'Regularized Coulomb friction law 2 parameters:\n'

        s += '   Regularized Coulomb friction law reads:\n'
        s += '                       C N |u|^(1/m)         \n'
        s += '      tau_b = -  ____________________________\n'
        s += '                   (|u| + (K*N)^m)^(1/m)     \n'
        s += '\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'C', 'friction coefficient [SI]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'm', 'm exponent'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'K', '(K * N) ^ m to be velocity controlling plastic limit'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.regcoulomb2 Class'
        return s

## ------------------------------------------------------
## friction.schoof
## ------------------------------------------------------
@class_registry.register_class
class schoof(class_registry.manage_state):
    '''
    friction.schoof Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.C = np.nan
        self.Cmax = np.nan
        self.m = np.nan
        self.coupling = 0
        self.effective_pressure = np.nan
        self.effective_pressure_limit = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        # See Brondex et al. 2019 https://www.the-cryosphere.net/13/177/2019/
        s = 'Schoof sliding law parameters:\n'

        s += '   Schoof\'s sliding law reads:\n'
        s += '                         C^2 |u_b|^(m-1)                \n'
        s += '      tau_b = - _____________________________   u_b   \n'
        s += '               (1+(C^2/(Cmax N))^1/m |u_b| )^m          \n'
        s += '\n'
        s += "{}\n".format(build_utils.fielddisplay(self, 'C', 'friction coefficient [SI]'))
        s += "{}\n".format(build_utils.fielddisplay(self, 'Cmax', 'Iken\'s bound (typically between 0.17 and 0.84) [SI]'))
        s += "{}\n".format(build_utils.fielddisplay(self, 'm', 'm exponent (generally taken as m = 1/n = 1/3)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'coupling', 'Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: used coupled model (not implemented yet)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        s += "{}\n".format(build_utils.fielddisplay(self, 'effective_pressure_limit', 'fNeff do not allow to fall below a certain limit: effective_pressure_limit*rho_ice*g*thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.schoof Class'
        return s

## ------------------------------------------------------
## friction.shakti
## ------------------------------------------------------
@class_registry.register_class
class shakti(class_registry.manage_state):
    '''
    friction.shakti Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.coefficient = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Basal shear stress parameters: Sigma_b = coefficient^2 * Neff * u_b\n'

        s += '(effective stress Neff = rho_ice * g * thickness + rho_water * g * (head - b))\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'coefficient', 'friction coefficient [SI]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.shakti Class'
        return s

## ------------------------------------------------------
## friction.waterlayer
## ------------------------------------------------------
@class_registry.register_class
class waterlayer(class_registry.manage_state):
    '''
    friction.waterlayer Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.coefficient = np.nan
        self.f = np.nan
        self.p = np.nan
        self.q = np.nan
        self.water_layer = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Basal shear stress parameters: tau_b = coefficient^2 * Neff ^r * |u_b|^(s - 1) * u_b * 1 / f(T)\n(effective stress Neff = rho_ice * g * thickness + rho_water * g * (bed + water_layer), r = q / p and s = 1 / p)\n'

        s = "{}\n".format(build_utils.fielddisplay(self, 'coefficient', 'frictiontemp coefficient [SI]'))
        s = "{}\n".format(build_utils.fielddisplay(self, 'f', 'f variable for effective pressure'))
        s = "{}\n".format(build_utils.fielddisplay(self, 'p', 'p exponent'))
        s = "{}\n".format(build_utils.fielddisplay(self, 'q', 'q exponent'))
        s = "{}\n".format(build_utils.fielddisplay(self, 'water_layer', 'water thickness at the base of the ice (m)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.waterlayer Class'
        return s

## ------------------------------------------------------
## friction.weertman
## ------------------------------------------------------
@class_registry.register_class
class weertman(class_registry.manage_state):
    '''
    friction.weertman Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.C = np.nan
        self.m = np.nan
        self.linearize = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Weertman sliding law parameters: Sigma_b = C^(- 1 / m) * |u_b|^(1 / m - 1) * u_b\n'

        s = "{}\n".format(build_utils.fielddisplay(self, 'C', 'friction coefficient [SI]'))
        s = "{}\n".format(build_utils.fielddisplay(self, 'm', 'm exponent'))
        s = "{}\n".format(build_utils.fielddisplay(self, 'linearize', '0: not linearized, 1: interpolated linearly, 2: constant per element (default is 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.weertman Class'
        return s