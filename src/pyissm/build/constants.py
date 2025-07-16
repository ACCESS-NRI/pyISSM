from . import build_utils
from . import class_registry

@class_registry.register_class
class constants(class_registry.manage_state):
    '''
    constants Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.g                      = 9.81
        self.omega                  = 7.292 * 1e-5
        self.yts                    = 365.0 * 24.0 * 3600.0
        self.referencetemperature   = 223.15
        self.gravitational_constant = 6.67259e-11

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   constants parameters:\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'g', 'gravitational acceleration [m/s^2]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'omega', 'angular velocity of Earth [rad/s]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'yts', 'number of seconds in a year [s/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'referencetemperature', 'reference temperature used in the enthalpy model [K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'gravitational_constant', 'Newtonian constant of gravitation [m^3/kg/s^2]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - constants Class'
        return s

