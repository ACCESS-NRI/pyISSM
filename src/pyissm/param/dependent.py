import numpy as np
from . import param_utils
from . import class_registry

@class_registry.register_class
class dependent(class_registry.manage_state):
    """
    Dependent variable parameters class for ISSM.

    This class encapsulates parameters for dependent variables in the ISSM (Ice Sheet System Model) framework.
    Dependent variables are outputs or responses that depend on independent variables and are typically
    used as objective functions in inverse problems or as outputs for sensitivity analysis.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    name : str, default=''
        Variable name (must match corresponding String).
    fos_reverse_index : float, default=nan
        Index for fos_reverse driver of ADOLC.
    exp : str, default=''
        File needed to compute dependent variable.
    segments : str, default='List of segments'
        Mass flux segments.
    index : int, default=-1
        Index parameter.
    nods : int, default=0
        Size parameter.

    Methods
    -------
    __init__(self, other=None)
        Initializes the dependent parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the dependent parameters.
    __str__(self)
        Returns a short string identifying the class.

    Notes
    -----
    This functionality is not yet fully implemented in the current version.

    Examples
    --------
    md.dependent = pyissm.param.dependent()
    md.dependent.name = 'Vel'
    md.dependent.exp = 'velocity_observations.exp'
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.name = ''
        self.fos_reverse_index = np.nan
        self.exp = ''
        self.segments = 'List of segments'
        self.index = -1
        self.nods = 0

        # Inherit matching fields from provided class
        super().__init__(other)

        # TODO: Implement check and adjustment for mass flux variable

    # Define repr
    def __repr__(self):
        s = '---------------------------------------\n'
        s += '****      NOT YET IMPLEMENTED      ****\n'
        s += '---------------------------------------\n\n'
        s += '   dependent variable:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'name', 'variable name (must match corresponding String)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'fos_reverse_index', 'index for fos_reverse driver of ADOLC'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'exp', 'file needed to compute dependent variable'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'segments', 'mass flux segments'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'index', '...'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'nods', '...'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - dependent Class'
        return s

