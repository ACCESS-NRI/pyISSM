import numpy as np
from . import param_utils
from . import class_registry

@class_registry.register_class
class independent(class_registry.manage_state):
    """
    Independent variable parameters class for ISSM.

    This class encapsulates parameters for independent variables in the ISSM (Ice Sheet System Model) framework.
    Independent variables are parameters that can be optimized or varied during inverse problems, 
    sensitivity analysis, or uncertainty quantification studies.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    name : str, default=''
        Variable name (must match corresponding String).
    type : str, default=''
        Type of variable ('vertex' or 'scalar').
    fos_forward_index : float, default=nan
        Index for fos_forward driver of ADOLC.
    fov_forward_indices : array, default=[]
        Indices for fov_forward driver of ADOLC.
    nods : int, default=0
        Size of independent variables.
    min_parameters : float, default=nan
        Absolute minimum acceptable value of the inversed parameter on each vertex.
    max_parameters : float, default=nan
        Absolute maximum acceptable value of the inversed parameter on each vertex.
    control_scaling_factor : float, default=nan
        Order of magnitude of each control (useful for multi-parameter optimization).
    control_size : int, default=1
        Number of timesteps.

    Methods
    -------
    __init__(self, other=None)
        Initializes the independent parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the independent parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.independent = pyissm.param.independent()
    md.independent.name = 'FrictionCoefficient'
    md.independent.type = 'vertex'
    md.independent.min_parameters = 1e-3
    md.independent.max_parameters = 100
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.name = ''
        self.type = ''
        self.fos_forward_index = np.nan
        self.fov_forward_indices = np.array([])
        self.nods = 0
        self.min_parameters = np.nan
        self.max_parameters = np.nan
        self.control_scaling_factor = np.nan
        self.control_size = 1

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   independent variable:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'name', 'variable name (must match corresponding String)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'type', 'type of variable (\'vertex\' or \'scalar\')'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'nods', 'size of independent variables'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'control_size', 'number of timesteps'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'min_parameters', 'absolute minimum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'max_parameters', 'absolute maximum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'control_scaling_factor', 'order of magnitude of each control (useful for multi-parameter optimization)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'fos_forward_index', 'index for fos_foward driver of ADOLC'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'fov_forward_indices', 'indices for fov_foward driver of ADOLC'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - independent Class'
        return s

