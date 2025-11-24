import collections
from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry

@class_registry.register_class
class private(class_registry.manage_state):
    """
    Private parameters class for ISSM.

    This class encapsulates internal parameters used by the ISSM (Ice Sheet System Model) framework.
    These parameters are not intended to be modified directly by users and are primarily used for
    managing model consistency, runtime information, mesh properties, and solution type.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    isconsistent : bool, default=1
        Indicates whether the model is self-consistent.
    runtimename : str, default=''
        Name of the run launched.
    bamg : collections.OrderedDict, default=OrderedDict()
        Structure with mesh properties constructed if BAMG is used to mesh the domain.
    solution : str, default=''
        Type of solution launched.

    Methods
    -------
    __init__(self, other=None)
        Initializes the private parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the private parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.private = pyissm.model.classes.private()
    md.private.runtimename = 'experiment_001'
    md.private.solution = 'StressBalance'
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.isconsistent = 1
        self.runtimename = ''
        self.bamg = collections.OrderedDict()
        self.solution = ''

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   private parameters -- do not change:\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'isconsistent', 'is model self consistent?'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'runtimename', 'name of the run launched'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'bamg', 'structure with mesh properties constructed if bamg is used to mesh the domain'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'solution', 'type of solution launched'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - private Class'
        return s

    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        # No checks
        return md
