from . import build_utils
from . import class_registry

@class_registry.register_class
class steadystate(class_registry.manage_state):
    """
    Steady state solution parameters class for ISSM.

    This class encapsulates parameters for configuring steady state simulations in the ISSM (Ice Sheet System Model) framework.
    It allows users to control convergence criteria and iteration limits for finding steady state solutions
    where the ice sheet geometry and flow fields reach equilibrium.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    reltol : float, default=0.01
        Relative tolerance criterion for convergence.
    maxiter : int, default=100
        Maximum number of iterations allowed.
    requested_outputs : str, default='List of requested outputs'
        Additional requested outputs for the steady state solution.

    Methods
    -------
    __init__(self, other=None)
        Initializes the steadystate parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the steadystate parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.steadystate = pyissm.build.steadystate()
    md.steadystate.reltol = 0.001
    md.steadystate.maxiter = 200
    md.steadystate.requested_outputs = ['IceVolume', 'IceVolumeAboveFloatation']
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.reltol = 0.01
        self.maxiter = 100
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   steadystate solution parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'reltol', 'relative tolerance criterion'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'maxiter', 'maximum number of iterations'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional requested outputs'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - steadystate Class'
        return s

