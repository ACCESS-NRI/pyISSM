import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class esa(class_registry.manage_state):
    '''
    esa Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.deltathickness = np.nan
        self.love_h = 0.
        self.love_l = 0.
        self.hemisphere = 0.
        self.degacc = 0.01
        self.requested_outputs = 'List of requested outputs'
        self.transitions = 'List of transitions'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   esa parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'deltathickness', 'thickness change: ice height equivalent [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'love_h', 'load Love number for radial displacement'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'love_l', 'load Love number for horizontal displaements'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'hemisphere', 'North-south, East-west components of 2-D horiz displacement vector:-1 south, 1 north'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'degacc', 'accuracy (default .01 deg) for numerical discretization of the Green''s functions'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'transitions', 'indices into parts of the mesh that will be icecaps'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested (default: EsaUmotion)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - esa Class'
        return s

