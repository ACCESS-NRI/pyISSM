import numpy as np
from . import param_utils
from . import class_registry

@class_registry.register_class
class sampling(class_registry.manage_state):
    """
    Sampling parameters class for ISSM.

    This class encapsulates parameters for stochastic sampling in the ISSM (Ice Sheet System Model) framework.
    It configures parameters for generating random fields using PDE-based operators and autoregressive processes,
    useful for uncertainty quantification and stochastic forcing applications.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    kappa : float, default=nan
        Coefficient of the identity operator in PDE operator (kappa^2 I - Laplacian)^(alpha/2)(tau).
    tau : float, default=0
        Scaling coefficient of the solution.
    beta : float, default=nan
        Coefficient in Robin boundary conditions (to be defined for robin = 1).
    phi : float, default=nan
        Temporal correlation factor for first-order autoregressive process X_t = phi X_{t-1} + noise (|phi|<1 for stationary process, phi = 1 for random walk process).
    alpha : float, default=2
        Exponent in PDE operator (default: 2.0, BiLaplacian covariance operator).
    robin : int, default=0
        Apply Robin boundary conditions (1 if applied and 0 for homogeneous Neumann boundary conditions).
    seed : int, default=-1
        Seed for pseudorandom number generator (given seed if >=0 and random seed if <0).
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested (not implemented yet).

    Methods
    -------
    __init__(self, other=None)
        Initializes the sampling parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the sampling parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.sampling = pyissm.build.sampling()
    md.sampling.kappa = 0.1
    md.sampling.alpha = 2.0
    md.sampling.phi = 0.9
    md.sampling.robin = 1
    md.sampling.beta = 0.5
    """

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
        s = '   Sampling parameters:\n'

        s += '      Parameters of PDE operator (kappa^2 I-Laplacian)^(alpha/2)(tau):\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'kappa', 'coefficient of the identity operator'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'tau', 'scaling coefficient of the solution (default: 1.0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'alpha', 'exponent in PDE operator, (default: 2.0, BiLaplacian covariance operator)'))
        s += '\n'
        s += '      Parameters of Robin boundary conditions nabla () \cdot normvec + beta ():\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'robin', 'Apply Robin boundary conditions (1 if applied and 0 for homogenous Neumann boundary conditions) (default: 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'beta', 'Coefficient in Robin boundary conditions (to be defined for robin = 1)'))
        s += '\n'
        s += '      Parameters for first-order autoregressive process (X_t = phi X_{t-1} + noise) (if transient):\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'phi', 'Temporal correlation factor (|phi|<1 for stationary process, phi = 1 for random walk process) (default 0)'))
        s += '\n'
        s += '      Other parameters of stochastic sampler:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'seed', 'Seed for pseudorandom number generator (given seed if >=0 and random seed if <0) (default: -1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested (not implemented yet)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - sampling Class'
        return s

