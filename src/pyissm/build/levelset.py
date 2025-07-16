import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class levelset(class_registry.manage_state):
    '''
    levelset Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.stabilization = 1
        self.spclevelset = np.nan
        self.reinit_frequency = 10
        self.kill_icebergs = 1
        self.migration_max = 1e12
        self.fe = 'P1'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Level-set parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'stabilization', '0: No Stabilization - No stabilization techniques applied.'))
        s += '{}\n'.format('                             1: Artificial Diffusivity - Most stable, but least accurate.')
        s += '{}\n'.format('                             2: Streamline Upwinding')
        s += '{}\n'.format('                             5: SUPG - Most accurate, but may be unstable in some applications.')
        s += '{}\n'.format(build_utils.fielddisplay(self, 'spclevelset', 'Levelset constraints (NaN means no constraint)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'reinit_frequency', 'Amount of time steps after which the levelset function in re-initialized'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'kill_icebergs', 'remove floating icebergs to prevent rigid body motions (1: true, 0: false)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'migration_max', 'maximum allowed migration rate (m/a)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'fe', 'Finite Element type: \'P1\' (default), or \'P2\''))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - levelset Class'
        return s

