import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class debris(class_registry.manage_state):
    '''
    debris Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.spcthickness = np.nan
        self.min_thickness = 0.
        self.stabilization = 2
        self.packingfraction = 0.01
        self.removalmodel = 0
        self.displacementmodel = 0
        self.max_displacementvelocity = 0.
        self.removal_slope_threshold = 0.
        self.removal_stress_threshold = 0.
        self.vertex_pairing = np.nan
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   debris solution parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self,'spcthickness','debris thickness constraints (NaN means no constraint) [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'min_thickness','minimum debris thickness allowed [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'packingfraction','fraction of debris covered in the ice'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'stabilization','0: no stabilization, 1: artificial diffusion, 2: streamline upwinding, 3: streamline upwind Petrov-Galerkin (SUPG)'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'removalmodel','frontal removal of debris. 0: no removal, 1: Slope-triggered debris removal, 2: driving-stress triggered debris removal'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'displacementmodel','debris displacement. 0: no displacement, 1: ...'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'max_displacementvelocity','maximum velocity of debris transport (v_ice + v_displacement) (m/a)'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'removal_slope_threshold','critical slope (degrees) for removalmodel (1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'removal_stress_threshold','critical stress (Pa) for removalmodel (2)'))

        s += '\n      {}\n'.format('Penalty options:')
        s += '{}\n'.format(build_utils.fielddisplay(self,'vertex_pairing','pairs of vertices that are penalized'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'requested_outputs','additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - debris Class'
        return s

