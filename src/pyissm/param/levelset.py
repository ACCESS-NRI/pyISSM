import numpy as np
from . import param_utils
from . import class_registry

@class_registry.register_class
class levelset(class_registry.manage_state):
    """
    Level-set method parameters class for ISSM.

    This class encapsulates parameters for configuring level-set method simulations in the ISSM (Ice Sheet System Model) framework.
    The level-set method is used to track moving boundaries such as ice fronts and calving fronts, 
    allowing for dynamic changes in ice sheet geometry during simulations.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    stabilization : int, default=1
        Stabilization method: 0=No Stabilization, 1=Artificial Diffusivity (most stable, least accurate), 2=Streamline Upwinding, 5=SUPG (most accurate, may be unstable).
    spclevelset : ndarray, default=nan
        Levelset constraints (NaN means no constraint).
    reinit_frequency : int, default=10
        Amount of time steps after which the levelset function is re-initialized.
    kill_icebergs : int, default=1
        Remove floating icebergs to prevent rigid body motions (1=true, 0=false).
    migration_max : float, default=1e12
        Maximum allowed migration rate [m/a].
    fe : str, default='P1'
        Finite Element type: 'P1' (default) or 'P2'.

    Methods
    -------
    __init__(self, other=None)
        Initializes the levelset parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the levelset parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.levelset = pyissm.param.levelset()
    md.levelset.stabilization = 5
    md.levelset.reinit_frequency = 5
    md.levelset.kill_icebergs = 0
    md.levelset.migration_max = 1e10
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.stabilization = 1
        self.spclevelset = np.nan
        self.reinit_frequency = 10
        self.kill_icebergs = 1
        self.migration_max = 1e12
        self.fe = 'P1'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Level-set parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'stabilization', '0: No Stabilization - No stabilization techniques applied.'))
        s += '{}\n'.format('                             1: Artificial Diffusivity - Most stable, but least accurate.')
        s += '{}\n'.format('                             2: Streamline Upwinding')
        s += '{}\n'.format('                             5: SUPG - Most accurate, but may be unstable in some applications.')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'spclevelset', 'Levelset constraints (NaN means no constraint)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'reinit_frequency', 'Amount of time steps after which the levelset function in re-initialized'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'kill_icebergs', 'remove floating icebergs to prevent rigid body motions (1: true, 0: false)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'migration_max', 'maximum allowed migration rate (m/a)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'fe', 'Finite Element type: \'P1\' (default), or \'P2\''))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - levelset Class'
        return s

