import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class mask(class_registry.manage_state):
    '''
    mask Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.ice_levelset = np.nan
        self.ocean_levelset = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   mask parameters:\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'ice_levelset', 'presence of ice if < 0, icefront position if = 0, no ice if > 0'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'ocean_levelset', 'presence of ocean if < 0, coastline/grounding line if = 0, no ocean if > 0'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - mask Class'
        return s

