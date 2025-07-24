import numpy as np
from . import param_utils
from . import class_registry

@class_registry.register_class
class thermal(class_registry.manage_state):
    """
    Thermal solution parameters class for ISSM.

    This class encapsulates parameters for configuring thermal simulations in the ISSM (Ice Sheet System Model) framework.
    It allows users to configure temperature constraints, stabilization methods, convergence criteria, and enthalpy formulations
    for solving the thermal evolution of ice sheets including temperate ice behavior.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    spctemperature : ndarray, default=nan
        Temperature constraints (NaN means no constraint) [K].
    penalty_threshold : float, default=0
        Threshold to declare convergence of thermal solution.
    stabilization : int, default=1
        Stabilization method: 0=no, 1=artificial_diffusivity, 2=SUPG.
    reltol : float, default=0.01
        Relative tolerance criterion for convergence.
    maxiter : int, default=100
        Maximum number of non-linear iterations.
    penalty_lock : int, default=0
        Stabilize unstable thermal constraints that keep zigzagging after n iterations (0=no stabilization).
    penalty_factor : int, default=3
        Penalty factor for constraint stabilization.
    isenthalpy : int, default=0
        Use an enthalpy formulation to include temperate ice.
    isdynamicbasalspc : int, default=0
        Enable dynamic setting of basal forcing. Required for enthalpy formulation.
    isdrainicecolumn : int, default=1
        Whether waterfraction drainage is enabled for enthalpy formulation.
    watercolumn_upperlimit : float, default=1000
        Upper limit of basal watercolumn for enthalpy formulation [m].
    fe : str, default='P1'
        Finite Element type: 'P1' (default), 'P1xP2'.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the thermal parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the thermal parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.thermal = pyissm.param.thermal()
    md.thermal.isenthalpy = 1
    md.thermal.stabilization = 2
    md.thermal.maxiter = 200
    md.thermal.reltol = 0.001
    """

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
        s += '{}\n'.format(param_utils.fielddisplay(self, 'spctemperature', 'temperature constraints (NaN means no constraint) [K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'stabilization', '0: no, 1: artificial_diffusivity, 2: SUPG'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'maxiter', 'maximum number of non linear iterations'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'reltol', 'relative tolerance criterion'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'penalty_lock', 'stabilize unstable thermal constraints that keep zigzagging after n iteration (default is 0, no stabilization)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'penalty_threshold', 'threshold to declare convergence of thermal solution (default is 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isenthalpy', 'use an enthalpy formulation to include temperate ice (default is 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isdynamicbasalspc', 'enable dynamic setting of basal forcing. required for enthalpy formulation (default is 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isdrainicecolumn', 'wether waterfraction drainage is enabled for enthalpy formulation (default is 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'watercolumn_upperlimit', 'upper limit of basal watercolumn for enthalpy formulation (default is 1000m)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'fe', 'Finite Element type: ''P1'' (default), ''P1xP2'''))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - thermal Class'
        return s

