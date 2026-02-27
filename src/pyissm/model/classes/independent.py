import numpy as np
from pyissm.model.classes import class_utils, class_registry

@class_registry.register_class
class independent(class_registry.manage_state):
    """
    Independent variable class for ISSM.

    This class contains parameters for independent variables in the ISSM framework.
    Independent variables are parameters that can be optimized or varied during inverse problems, 
    sensitivity analysis, or uncertainty quantification studies.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    name : :class:`str`, default=''
        Variable name (must match corresponding String).
    type : :class:`str`, default=''
        Type of variable ('vertex' or 'scalar').
    fos_forward_index : :class:`float`, default=np.nan
        Index for fos_forward driver of ADOLC.
    fov_forward_indices : :class:`numpy.ndarray`, default=np.array([])
        Indices for fov_forward driver of ADOLC.
    nods : :class:`int`, default=0
        Size of independent variables.
    min_parameters : :class:`float`, default=np.nan
        Absolute minimum acceptable value of the inversed parameter on each vertex.
    max_parameters : :class:`float`, default=np.nan
        Absolute maximum acceptable value of the inversed parameter on each vertex.
    control_scaling_factor : :class:`float`, default=1.0
        Order of magnitude of each control (useful for multi-parameter optimization).
    control_size : :class:`int`, default=1
        Number of timesteps.

    Examples
    --------
    .. code-block:: python
    
        >>> md.independent = pyissm.model.classes.independent()
        >>> md.independent.name = 'FrictionCoefficient'
        >>> md.independent.type = 'vertex'
        >>> md.independent.min_parameters = 1e-3
        >>> md.independent.max_parameters = 100
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
        self.control_scaling_factor = 1.0
        self.control_size = 1

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   independent variable:\n'

        s += '{}\n'.format(class_utils._field_display(self, 'name', 'variable name (must match corresponding String)'))
        s += '{}\n'.format(class_utils._field_display(self, 'type', 'type of variable (\'vertex\' or \'scalar\')'))
        s += '{}\n'.format(class_utils._field_display(self, 'nods', 'size of independent variables'))
        s += '{}\n'.format(class_utils._field_display(self, 'control_size', 'number of timesteps'))
        s += '{}\n'.format(class_utils._field_display(self, 'min_parameters', 'absolute minimum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(class_utils._field_display(self, 'max_parameters', 'absolute maximum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(class_utils._field_display(self, 'control_scaling_factor', 'order of magnitude of each control (useful for multi-parameter optimization)'))
        s += '{}\n'.format(class_utils._field_display(self, 'fos_forward_index', 'index for fos_foward driver of ADOLC'))
        s += '{}\n'.format(class_utils._field_display(self, 'fov_forward_indices', 'indices for fov_foward driver of ADOLC'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - independent Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [independent] parameters.

        Parameters
        ----------
        md : :class:`pyissm.model.Model`
            The model object to check.
        solution : :class:`pyissm.model.solution`
            The solution object to check.
        analyses : list of :class:`str`
            List of analyses to check consistency for.

        Returns 
        -------
        md : :class:`pyissm.model.Model`
            The model object with any consistency errors noted.
        """

        if not np.isnan(self.fos_forward_index):
            if self.nods == 0:
                raise TypeError('pyissm.model.classes.independent.check_consistency: nods should be set to the size of the independent variable')

        if len(self.fov_forward_indices) > 0:
            if self.nods == 0:
                raise TypeError('pyissm.model.classes.independent.check_consistency: nods should be set to the size of the independent variable')
            
            class_utils.check_field(md, field = self.fov_forward_indices, ge = 1, le = self.nods, message = "pyissm.model.classes.independent.check_consistency: fov_forward_indices should be between 1 and nods (inclusive).")
        
        md = class_utils.check_field(md, field = self.control_scaling_factor, scalar = True, gt = 0., allow_nan = False, allow_inf = False, message = "pyissm.model.classes.independent.check_consistency: control_scaling_factor should be a positive scalar value.")

        return md