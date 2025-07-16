import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class damage(class_registry.manage_state):
    '''
    damage Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.isdamage = 0
        self.D = np.nan
        self.law = np.nan
        self.spcdamage = np.nan
        self.max_damage = np.nan
        self.stabilization = np.nan
        self.maxiter = np.nan
        self.elementinterp = ''
        self.stress_threshold = np.nan
        self.stress_ubound = np.nan
        self.kappa = np.nan
        self.c1 = np.nan
        self.c2 = np.nan
        self.c3 = np.nan
        self.c4 = np.nan
        self.healing = np.nan
        self.equiv_stress = np.nan
        self.requested_outputs = 'List requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Damage:\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isdamage', 'is damage mechanics being used? [0 (default) or 1]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, "D", "damage tensor (scalar for now)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "law", "damage law ['0: analytical', '1: pralong']"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "spcdamage", "damage constraints (NaN means no constraint)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "max_damage", "maximum possible damage (0 <=max_damage < 1)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "stabilization", "0: no stabilization, 1: artificial diffusion, 2: SUPG (not working), 4: flux corrected transport"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "maxiter", "maximum number of non linear iterations"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "elementinterp", "interpolation scheme for finite elements [''P1'', ''P2'']"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "stress_threshold", "stress threshold for damage initiation (Pa)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "stress_ubound", "stress upper bound for damage healing (Pa)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "kappa", "ductility parameter for stress softening and damage [ > 1]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "c1", "damage parameter 1 "))
        s += '{}\n'.format(build_utils.fielddisplay(self, "c2", "damage parameter 2 "))
        s += '{}\n'.format(build_utils.fielddisplay(self, "c3", "damage parameter 3 "))
        s += '{}\n'.format(build_utils.fielddisplay(self, "c4", "damage parameter 4 "))
        s += '{}\n'.format(build_utils.fielddisplay(self, "healing", "damage healing parameter"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "equiv_stress", "0: von Mises, 1: max principal"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - damage Class'
        return s

