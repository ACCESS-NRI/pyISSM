import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class radaroverlay(class_registry.manage_state):
    '''
    radaroverlay Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.pwr = np.nan
        self.x = np.nan
        self.y = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   radaroverlay parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'pwr', 'radar power image (matrix)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'x', 'corresponding x coordinates [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'y', 'corresponding y coordinates [m]'))

        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - radaroverlay Class'
        return s

