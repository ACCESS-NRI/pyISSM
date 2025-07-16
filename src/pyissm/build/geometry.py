import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class geometry(class_registry.manage_state):
    '''
    geometry Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.surface = np.nan
        self.thickness = np.nan
        self.base = np.nan
        self.bed = np.nan
        self.hydrostatic_ratio = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   geometry parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'surface', 'ice upper surface elevation [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thickness', 'ice thickness [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'base', 'ice base elevation [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'bed', 'bed elevation [m]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - geometry Class'
        return s

