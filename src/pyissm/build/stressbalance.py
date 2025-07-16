import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class stressbalance(class_registry.manage_state):
    '''
    stressbalance Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.spcvx = np.nan
        self.spcvy = np.nan
        self.spcvx_base = np.nan
        self.spcvy_base = np.nan
        self.spcvx_shear = np.nan
        self.spcvy_shear = np.nan
        self.spcvz = np.nan
        self.restol = pow(10, -4)
        self.reltol = 0.01
        self.abstol = 10
        self.ishydrologylayer = 0
        self.isnewton = 0
        self.FSreconditioning = pow(10, 13)
        #self.icefront = np.nan -- no longer in use
        self.maxiter = 100
        self.shelf_dampening = 0
        self.vertex_pairing = np.nan
        self.penalty_factor = 3
        self.rift_penalty_lock = 10
        self.rift_penalty_threshold = 0
        self.referential = np.nan
        self.loadingforce = np.nan
        self.requested_outputs = 'List of requested outputs' # Default = ['default'] (Vx, Vy, Vel, Pressure)

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   StressBalance solution parameters:\n'

        s += '      Convergence criteria:\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'restol', 'mechanical equilibrium residual convergence criterion'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'reltol', 'velocity relative convergence criterion, NaN: not applied'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'abstol', 'velocity absolute convergence criterion, NaN: not applied'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isnewton', '0: Picard\'s fixed point, 1: Newton\'s method, 2: hybrid'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'maxiter', 'maximum number of nonlinear iterations'))
        s += '\n'
        s += '      boundary conditions:\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'spcvx', 'x-axis velocity constraint (NaN means no constraint) [m / yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'spcvy', 'y-axis velocity constraint (NaN means no constraint) [m / yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'spcvz', 'z-axis velocity constraint (NaN means no constraint) [m / yr]'))
        s += '\n'
        s += '      MOLHO boundary conditions:\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'spcvx_base', 'x-axis basal velocity constraint (NaN means no constraint) [m / yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'spcvy_base', 'y-axis basal velocity constraint (NaN means no constraint) [m / yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'spcvx_shear', 'x-axis shear velocity constraint (NaN means no constraint) [m / yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'spcvy_shear', 'y-axis shear velocity constraint (NaN means no constraint) [m / yr]'))
        s += '\n'
        s += '      Rift options:\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rift_penalty_threshold', 'threshold for instability of mechanical constraints'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rift_penalty_lock', 'number of iterations before rift penalties are locked'))
        s += '\n'
        s += '      Penalty options:\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'penalty_factor', 'offset used by penalties: penalty = Kmax * 10^offset'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vertex_pairing', 'pairs of vertices that are penalized'))
        s += '\n'
        s += '      Hydrology layer:\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'ishydrologylayer', '(SSA only) 0: no subglacial hydrology layer in driving stress, 1: hydrology layer in driving stress'));
        s += '\n'
        s += '      Other:\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'shelf_dampening', 'use dampening for floating ice ? Only for FS model'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'FSreconditioning', 'multiplier for incompressibility equation. Only for FS model'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'referential', 'local referential'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'loadingforce', 'loading force applied on each point [N/m^3]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - stressbalance Class'
        return s

