import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class debris(class_registry.manage_state):
    """
    Debris transport parameters class for ISSM.

    This class encapsulates parameters for debris transport modeling in the ISSM (Ice Sheet System Model) framework.
    Debris transport simulates the movement and accumulation of rock debris on glacier surfaces,
    which affects surface albedo, melting rates, and overall ice dynamics.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    spcthickness : ndarray, default=nan
        Debris thickness constraints (NaN means no constraint) [m].
    min_thickness : float, default=0.0
        Minimum debris thickness allowed [m].
    stabilization : int, default=2
        Stabilization method: 0=no stabilization, 1=artificial diffusion, 2=streamline upwinding, 3=streamline upwind Petrov-Galerkin (SUPG).
    packingfraction : float, default=0.01
        Fraction of debris covered in the ice.
    removalmodel : int, default=0
        Frontal removal of debris: 0=no removal, 1=Slope-triggered debris removal, 2=driving-stress triggered debris removal.
    displacementmodel : int, default=0
        Debris displacement: 0=no displacement, 1=...
    max_displacementvelocity : float, default=0.0
        Maximum velocity of debris transport (v_ice + v_displacement) [m/a].
    removal_slope_threshold : float, default=0.0
        Critical slope [degrees] for removalmodel (1).
    removal_stress_threshold : float, default=0.0
        Critical stress [Pa] for removalmodel (2).
    vertex_pairing : float, default=nan
        Pairs of vertices that are penalized.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the debris parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the debris parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.debris = pyissm.build.debris()
    md.debris.min_thickness = 0.001
    md.debris.packingfraction = 0.02
    md.debris.stabilization = 2
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.spcthickness = np.nan
        self.min_thickness = 0.
        self.stabilization = 2
        self.packingfraction = 0.01
        self.removalmodel = 0
        self.displacementmodel = 0
        self.max_displacementvelocity = 0.
        self.removal_slope_threshold = 0.
        self.removal_stress_threshold = 0.
        self.vertex_pairing = np.nan
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   debris solution parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self,'spcthickness','debris thickness constraints (NaN means no constraint) [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'min_thickness','minimum debris thickness allowed [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'packingfraction','fraction of debris covered in the ice'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'stabilization','0: no stabilization, 1: artificial diffusion, 2: streamline upwinding, 3: streamline upwind Petrov-Galerkin (SUPG)'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'removalmodel','frontal removal of debris. 0: no removal, 1: Slope-triggered debris removal, 2: driving-stress triggered debris removal'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'displacementmodel','debris displacement. 0: no displacement, 1: ...'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'max_displacementvelocity','maximum velocity of debris transport (v_ice + v_displacement) (m/a)'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'removal_slope_threshold','critical slope (degrees) for removalmodel (1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'removal_stress_threshold','critical stress (Pa) for removalmodel (2)'))

        s += '\n      {}\n'.format('Penalty options:')
        s += '{}\n'.format(build_utils.fielddisplay(self,'vertex_pairing','pairs of vertices that are penalized'))
        s += '{}\n'.format(build_utils.fielddisplay(self,'requested_outputs','additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - debris Class'
        return s

