from . import param_utils
from . import class_registry

@class_registry.register_class
class constants(class_registry.manage_state):
    """
    Physical constants class for ISSM.

    This class encapsulates fundamental physical constants used in the ISSM (Ice Sheet System Model) framework.
    It provides standardized values for gravitational acceleration, Earth's rotation, time conversions, 
    and other physical constants required for ice sheet modeling calculations.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    g : float, default=9.81
        Gravitational acceleration [m/s^2].
    omega : float, default=7.292e-5
        Angular velocity of Earth [rad/s].
    yts : float, default=31536000.0
        Number of seconds in a year [s/yr] (365.0 * 24.0 * 3600.0).
    referencetemperature : float, default=223.15
        Reference temperature used in the enthalpy model [K].
    gravitational_constant : float, default=6.67259e-11
        Newtonian constant of gravitation [m^3/kg/s^2].

    Methods
    -------
    __init__(self, other=None)
        Initializes the physical constants, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the constants.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.constants = pyissm.build.constants()
    md.constants.g = 9.81
    md.constants.yts = 365.25 * 24.0 * 3600.0
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.g = 9.81
        self.omega = 7.292 * 1e-5
        self.yts = 365.0 * 24.0 * 3600.0
        self.referencetemperature = 223.15
        self.gravitational_constant = 6.67259e-11

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   constants parameters:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'g', 'gravitational acceleration [m/s^2]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'omega', 'angular velocity of Earth [rad/s]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'yts', 'number of seconds in a year [s/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'referencetemperature', 'reference temperature used in the enthalpy model [K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'gravitational_constant', 'Newtonian constant of gravitation [m^3/kg/s^2]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - constants Class'
        return s

