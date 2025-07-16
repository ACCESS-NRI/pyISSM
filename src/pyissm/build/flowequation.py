import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class flowequation(class_registry.manage_state):
    '''
    flowequation Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.isSIA = 0
        self.isL1L2 = 0
        self.isSSA = 0
        self.isMOLHO = 0
        self.isHO = 0
        self.isFS = 0
        self.isNitscheBC = 0
        self.FSNitscheGamma = 1e6
        self.fe_SSA = 'P1'
        self.fe_HO = 'P1'
        self.fe_FS = 'MINIcondensed'
        self.augmented_lagrangian_r = 1
        self.augmented_lagrangian_rhop = 1
        self.augmented_lagrangian_rlambda = 1
        self.augmented_lagrangian_rholambda = 1
        self.XTH_theta = 0
        self.vertex_equation = np.nan
        self.element_equation = np.nan
        self.borderSSA = np.nan
        self.borderHO = np.nan
        self.borderFS = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   flow equation parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'isSIA', "is the Shallow Ice Approximation (SIA) used?"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isSSA', "is the Shelfy-Stream Approximation (SSA) used?"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isL1L2', "are L1L2 equations used?"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isMOLHO', "are MOno-layer Higher-Order (MOLHO) equations used?"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isHO', "is the Higher-Order (HO) approximation used?"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isFS', "are the Full-FS (FS) equations used?"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isNitscheBC', "is weakly imposed condition used?"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'FSNitscheGamma', "Gamma value for the Nitsche term (default: 1e6)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'fe_SSA', "Finite Element for SSA: 'P1', 'P1bubble' 'P1bubblecondensed' 'P2'"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'fe_HO', "Finite Element for HO:  'P1', 'P1bubble', 'P1bubblecondensed', 'P1xP2', 'P2xP1', 'P2', 'P2bubble', 'P1xP3', 'P2xP4'"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'fe_FS', "Finite Element for FS:  'P1P1' (debugging only) 'P1P1GLS' 'MINIcondensed' 'MINI' 'TaylorHood' 'LATaylorHood' 'XTaylorHood'"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vertex_equation', "flow equation for each vertex"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'element_equation', "flow equation for each element"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'borderSSA', "vertices on SSA's border (for tiling)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'borderHO', "vertices on HO's border (for tiling)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'borderFS', "vertices on FS' border (for tiling)"))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - flowequation Class'
        return s

