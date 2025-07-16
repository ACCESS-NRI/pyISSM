import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class independent(class_registry.manage_state):
    '''
    independent Class definition
    '''

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

        s += '{}\n'.format(build_utils.fielddisplay(self, 'name', 'variable name (must match corresponding String)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'type', 'type of variable (\'vertex\' or \'scalar\')'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'nods', 'size of independent variables'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'control_size', 'number of timesteps'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'min_parameters', 'absolute minimum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'max_parameters', 'absolute maximum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'control_scaling_factor', 'order of magnitude of each control (useful for multi-parameter optimization)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'fos_forward_index', 'index for fos_foward driver of ADOLC'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'fov_forward_indices', 'indices for fov_foward driver of ADOLC'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - independent Class'
        return s

