import numpy as np
from . import build_utils
from . import class_registry

## ------------------------------------------------------
## calving.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    '''
    calving.default Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.calvingrate = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'calvingrate', 'calving rate at given location [m/a]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.default Class'
        return s

## ------------------------------------------------------
## calving.crevassedepth
## ------------------------------------------------------
@class_registry.register_class
class crevassedepth(class_registry.manage_state):
    '''
    calving.crevassedepth Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.crevasse_opening_stress = 1.
        self.crevasse_threshold = 1.
        self.water_height = 0.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving Pi parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'crevasse_opening_stress', '0: stress only in the ice-flow direction, 1: max principal'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'crevasse_threshold', 'ratio of full thickness to calve (e.g. 0.75 is for 75% of the total ice thickness)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'water_height', 'water height in the crevasse [m]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.crevassedepth Class'
        return s

## ------------------------------------------------------
## calving.dev
## ------------------------------------------------------
@class_registry.register_class
class dev(class_registry.manage_state):
    '''
    calving.dev Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.stress_threshold_groundedice = 1e6
        self.stress_threshold_floatingice = 150e3

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving Pi parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'stress_threshold_groundedice', 'sigma_max applied to grounded ice only [Pa]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'stress_threshold_floatingice', 'sigma_max applied to floating ice only [Pa]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.dev Class'
        return s

## ------------------------------------------------------
## calving.levermann
## ------------------------------------------------------
@class_registry.register_class
class levermann(class_registry.manage_state):
    '''
    calving.levermann Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.coeff = 2e13

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving Levermann parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'coeff', 'proportionality coefficient in Levermann model'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.levermann Class'
        return s

## ------------------------------------------------------
## calving.minthickness
## ------------------------------------------------------
@class_registry.register_class
class minthickness(class_registry.manage_state):
    '''
    calving.minthickness Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.min_thickness = 100.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving Minimum thickness:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'min_thickness', 'minimum thickness below which no ice is allowed'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.minthickness Class'
        return s

## ------------------------------------------------------
## calving.parameterization
## ------------------------------------------------------
@class_registry.register_class
class parameterization(class_registry.manage_state):
    '''
    calving.parameterization Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.min_thickness = 0.
        self.use_param = 0.
        self.theta = 0.
        self.alpha = 0.
        self.xoffset = 0.
        self.yoffset = 0.
        self.vel_upperbound = 6000.
        self.vel_threshold = 0.
        self.vel_lowerbound = 0.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving test parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'min_thickness', 'minimum thickness below which no ice is allowed [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'use_param', '-1 - just use frontal ablation rate, 0 - f(x) = y_{o} + \alpha (x+x_{o}), 1 - f(x)=y_{o}-\frac{\theta}{2}\tanh(\alpha(x+x_{o})), 2 - tanh(thickness), 3 - tanh(normalized vel), 4 - tanh(truncated vel), 5 - linear(truncated vel)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'theta', 'the amplifier'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'alpha', 'the slope'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'xoffset', 'offset in x-axis'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'yoffset', 'offset in y-axis'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vel_lowerbound', 'lowerbound of ice velocity to reduce the calving rate [m/a]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vel_threshold', 'threshold of ice velocity to reduce the calving rate [m/a]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vel_upperbound', 'upperbound of ice velocity to reduce the calving rate [m/a]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.parameterization Class'
        return s

## ------------------------------------------------------
## calving.vonmises
## ------------------------------------------------------
@class_registry.register_class
class vonmises(class_registry.manage_state):
    '''
    calving.vonmises Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.stress_threshold_groundedice = 0
        self.stress_threshold_floatingice = 0
        self.min_thickness = 0.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving VonMises parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'stress_threshold_groundedice', 'sigma_max applied to grounded ice only [Pa]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'stress_threshold_floatingice', 'sigma_max applied to floating ice only [Pa]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'min_thickness', 'minimum thickness below which no ice is allowed [m]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.vonmises Class'
        return s