from . import build_utils
from . import class_registry

@class_registry.register_class
class groundingline(class_registry.manage_state):
    """
    Grounding line migration parameters class for ISSM.

    This class encapsulates parameters for configuring grounding line migration in the ISSM (Ice Sheet System Model) framework.
    It controls how the grounding line (boundary between grounded and floating ice) moves during simulations,
    including migration methods and interpolation schemes for friction and melting on partially floating elements.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    migration : str, default='SubelementMigration'
        Type of grounding line migration: 'SoftMigration', 'SubelementMigration', 'AggressiveMigration', 'Contact', 'None'.
    friction_interpolation : str, default='SubelementFriction1'
        Type of friction interpolation on partially floating elements: 'SubelementFriction1', 'SubelementFriction2', 'NoFrictionOnPartiallyFloating'.
    melt_interpolation : str, default='NoMeltOnPartiallyFloating'
        Type of melt interpolation on partially floating elements: 'SubelementMelt1', 'SubelementMelt2', 'IntrusionMelt', 'NoMeltOnPartiallyFloating', 'FullMeltOnPartiallyFloating'.
    intrusion_distance : float, default=0
        Distance of seawater intrusion from grounding line [m].
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested for grounding line analysis.

    Methods
    -------
    __init__(self, other=None)
        Initializes the groundingline parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the groundingline parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.groundingline = pyissm.build.groundingline()
    md.groundingline.migration = 'AggressiveMigration'
    md.groundingline.friction_interpolation = 'SubelementFriction2'
    md.groundingline.melt_interpolation = 'SubelementMelt1'
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.migration = 'SubelementMigration'
        self.friction_interpolation = 'SubelementFriction1'
        self.melt_interpolation = 'NoMeltOnPartiallyFloating'
        self.intrusion_distance = 0
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   grounding line migration parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'migration', 'type of grounding line migration: \'SoftMigration\', \'SubelementMigration\', \'AggressiveMigration\', \'Contact\', \'None\''))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'friction_interpolation', 'type of friction interpolation on partially floating elements: ''SubelementFriction1'', ''SubelementFriction2'', ''NoFrictionOnPartiallyFloating'''))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'melt_interpolation', 'type of melt interpolation on partially floating elements: \'SubelementMelt1\', \'SubelementMelt2\', \'IntrusionMelt\', \'NoMeltOnPartiallyFloating\', \'FullMeltOnPartiallyFloating\''))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'intrusion_distance', 'distance of seawater intrusion from grounding line [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - groundingline Class'
        return s

