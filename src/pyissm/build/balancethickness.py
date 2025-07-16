import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class balancethickness(class_registry.manage_state):
    '''
    balancethickness Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.spcthickness = np.nan
        self.thickening_rate = np.nan
        self.stabilization = 0
        self.omega = np.nan
        self.slopex = np.nan
        self.slopey = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   balance thickness solution parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'spcthickness', 'thickness constraints (NaN means no constraint) [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thickening_rate', 'ice thickening rate used in the mass conservation (dh / dt) [m / yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'stabilization', "0: None, 1: SU, 2: SSA's artificial diffusivity, 3:DG"))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - balancethickness Class'
        return s

