import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class surfaceload(class_registry.manage_state):
    '''
    surfaceload Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.icethicknesschange = np.nan
        self.waterheightchange = np.nan
        self.other = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surfaceload:\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'icethicknesschange', 'thickness change: ice height equivalent [mIce/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'waterheightchange', 'water height change: water height equivalent [mWater/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'other', 'other loads (sediments) [kg/m^2/yr]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - surfaceload Class'
        return s

