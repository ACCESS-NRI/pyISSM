import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class offlinesolidearthsolution(class_registry.manage_state):
    '''
    offlinesolidearthsolution Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.displacementeast = np.nan
        self.displacementnorth = np.nan
        self.displacementup = np.nan
        self.geoid = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    def __repr__(self):
        s = '         units for time series is (yr)\n       external: offlinesolidearth solution\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'displacementeast', 'solid-Earth Eastwards bedrock displacement time series (m)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'displacementnorth', 'solid-Earth Northwards bedrock displacement time series (m)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'displacementup', 'solid-Earth bedrock uplift time series (m)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'geoid', 'solid-Earth geoid time series (m)'))
        return s

    ## Define class string
    def __str__(self):
        s = 'ISSM - offlinesolidearthsolution Class'
        return s

