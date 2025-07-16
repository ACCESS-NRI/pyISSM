import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class thermal(class_registry.manage_state):
    '''
    thermal Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.spctemperature = np.nan
        self.penalty_threshold = 0
        self.stabilization = 1
        self.reltol = 0.01
        self.maxiter = 100
        self.penalty_lock = 0
        self.penalty_factor = 3
        self.isenthalpy = 0
        self.isdynamicbasalspc = 0
        self.isdrainicecolumn = 1
        self.watercolumn_upperlimit = 1000
        self.fe = 'P1'
        self.requested_outputs = 'List of requested outputs' # Default = ['default'] - Default varies depending on isenthalpy

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Thermal solution parameters:\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'spctemperature', 'temperature constraints (NaN means no constraint) [K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'stabilization', '0: no, 1: artificial_diffusivity, 2: SUPG'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'maxiter', 'maximum number of non linear iterations'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'reltol', 'relative tolerance criterion'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'penalty_lock', 'stabilize unstable thermal constraints that keep zigzagging after n iteration (default is 0, no stabilization)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'penalty_threshold', 'threshold to declare convergence of thermal solution (default is 0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isenthalpy', 'use an enthalpy formulation to include temperate ice (default is 0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isdynamicbasalspc', 'enable dynamic setting of basal forcing. required for enthalpy formulation (default is 0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isdrainicecolumn', 'wether waterfraction drainage is enabled for enthalpy formulation (default is 1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'watercolumn_upperlimit', 'upper limit of basal watercolumn for enthalpy formulation (default is 1000m)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'fe', 'Finite Element type: ''P1'' (default), ''P1xP2'''))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - thermal Class'
        return s

