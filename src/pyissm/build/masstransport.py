import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class masstransport(class_registry.manage_state):
    '''
    masstransport Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.spcthickness = np.nan
        self.isfreesurface = 0
        self.min_thickness = 1.
        self.hydrostatic_adjustment = 'Absolute'
        self.stabilization = 1
        self.vertex_pairing = np.nan
        self.penalty_factor = 3
        self.requested_outputs = 'List of requested outputs' # Default = ['default'] (Thickness, surface, base)

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Masstransport solution parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'spcthickness', 'thickness constraints (NaN means no constraint) [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isfreesurface', 'do we use free surfaces (FS only) or mass conservation'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'min_thickness', 'minimum ice thickness allowed [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'hydrostatic_adjustment', 'adjustment of ice shelves surface and bed elevations: ''Incremental'' or ''Absolute'' '))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'stabilization', '0: no stabilization, 1: artificial diffusion, 2: streamline upwinding, 3: discontinuous Galerkin, 4: flux corrected transport, 5: streamline upwind Petrov-Galerkin (SUPG)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - masstransport Class'
        return s

