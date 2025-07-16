import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class sampling(class_registry.manage_state):
    '''
    sampling Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.kappa = np.nan
        self.tau = 0
        self.beta = np.nan
        self.phi = np.nan
        self.alpha = 2
        self.robin = 0
        self.seed = -1
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Sampling parameters::\n'

        s += '      Parameters of PDE operator (kappa^2 I-Laplacian)^(alpha/2)(tau):\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'kappa', 'coefficient of the identity operator'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'tau', 'scaling coefficient of the solution (default: 1.0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'alpha', 'exponent in PDE operator, (default: 2.0, BiLaplacian covariance operator)'))
        s += '\n'
        s += '      Parameters of Robin boundary conditions nabla () \cdot normvec + beta ():\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'robin', 'Apply Robin boundary conditions (1 if applied and 0 for homogenous Neumann boundary conditions) (default: 0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'beta', 'Coefficient in Robin boundary conditions (to be defined for robin = 1)'))
        s += '\n'
        s += '      Parameters for first-order autoregressive process (X_t = phi X_{t-1} + noise) (if transient):\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'phi', 'Temporal correlation factor (|phi|<1 for stationary process, phi = 1 for random walk process) (default 0)'))
        s += '\n'
        s += '      Other parameters of stochastic sampler:\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'seed', 'Seed for pseudorandom number generator (given seed if >=0 and random seed if <0) (default: -1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested (not implemented yet)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - sampling Class'
        return s

