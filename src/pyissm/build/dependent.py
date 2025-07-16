import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class dependent(class_registry.manage_state):
    '''
    dependent Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.name = ''
        self.fos_reverse_index = np.nan
        self.exp = ''
        self.segments = 'List of segments'
        self.index = -1
        self.nods = 0

        # Inherit matching fields from provided class
        super().__init__(other)

        # TODO: Implement check and adjustment for mass flux variable

    # Define repr
    def __repr__(self):
        s = '---------------------------------------\n'
        s += '****      NOT YET IMPLEMENTED      ****\n'
        s += '---------------------------------------\n\n'
        s += '   dependent variable:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'name', 'variable name (must match corresponding String)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'fos_reverse_index', 'index for fos_reverse driver of ADOLC'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'exp', 'file needed to compute dependent variable'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'segments', 'mass flux segments'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'index', '...'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'nods', '...'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - dependent Class'
        return s

