import numpy as np
from . import param_utils
from . import class_registry

@class_registry.register_class
class masstransport(class_registry.manage_state):
    """
    Mass transport solution parameters class for ISSM.

    This class encapsulates parameters for configuring mass transport simulations in the ISSM (Ice Sheet System Model) framework.
    It controls ice thickness evolution, free surface behavior, stabilization methods, and hydrostatic adjustments
    for both grounded and floating ice.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    spcthickness : ndarray, default=nan
        Thickness constraints (NaN means no constraint) [m].
    isfreesurface : int, default=0
        Do we use free surfaces (FS only) or mass conservation.
    min_thickness : float, default=1.0
        Minimum ice thickness allowed [m].
    hydrostatic_adjustment : str, default='Absolute'
        Adjustment of ice shelves surface and bed elevations: 'Incremental' or 'Absolute'.
    stabilization : int, default=1
        Stabilization method: 0=no stabilization, 1=artificial diffusion, 2=streamline upwinding, 3=discontinuous Galerkin, 4=flux corrected transport, 5=streamline upwind Petrov-Galerkin (SUPG).
    vertex_pairing : float, default=nan
        Vertex pairing parameter. Used during consistency checks.
    penalty_factor : float, default=3
        Penalty factor for constraint enforcement.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the masstransport parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the masstransport parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.masstransport = pyissm.param.masstransport()
    md.masstransport.min_thickness = 10.0
    md.masstransport.stabilization = 2
    md.masstransport.hydrostatic_adjustment = 'Incremental'
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.spcthickness = np.nan
        self.isfreesurface = 0
        self.min_thickness = 1.
        self.hydrostatic_adjustment = 'Absolute'
        self.stabilization = 1
        self.vertex_pairing = np.nan
        self.penalty_factor = 3
        self.requested_outputs = 'List of requested outputs' # Default = ['default'] (Thickness, surface, base)

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Masstransport solution parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'spcthickness', 'thickness constraints (NaN means no constraint) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isfreesurface', 'do we use free surfaces (FS only) or mass conservation'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'min_thickness', 'minimum ice thickness allowed [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'hydrostatic_adjustment', 'adjustment of ice shelves surface and bed elevations: ''Incremental'' or ''Absolute'' '))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'stabilization', '0: no stabilization, 1: artificial diffusion, 2: streamline upwinding, 3: discontinuous Galerkin, 4: flux corrected transport, 5: streamline upwind Petrov-Galerkin (SUPG)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - masstransport Class'
        return s

