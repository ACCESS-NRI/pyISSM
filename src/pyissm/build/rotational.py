from . import build_utils
from . import class_registry

@class_registry.register_class
class rotational(class_registry.manage_state):
    """
    Rotational parameters class for ISSM.

    This class encapsulates parameters related to Earth's rotational properties for use in the ISSM (Ice Sheet System Model) framework.
    It stores moments of inertia and angular velocity values used in rotational calculations.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    equatorialmoi : float, default=8.0077e37
        Mean equatorial moment of inertia [kg m^2].
    polarmoi : float, default=8.0345e37
        Polar moment of inertia [kg m^2].
    angularvelocity : float, default=7.2921e-5
        Mean rotational velocity of earth [per second].

    Methods
    -------
    __init__(self, other=None)
        Initializes the rotational parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the rotational parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.rotational = pyissm.build.rotational()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.equatorialmoi = 8.0077 * pow(10, 37)
        self.polarmoi = 8.0345 * pow(10, 37)
        self.angularvelocity = 7.2921 * pow(10, -5)

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   rotational parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'equatorialmoi', 'mean equatorial moment of inertia [kg m^2]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'polarmoi', 'polar moment of inertia [kg m^2]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'angularvelocity', 'mean rotational velocity of earth [per second]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - rotational Class'
        return s

