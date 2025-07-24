from . import param_utils
from . import class_registry

@class_registry.register_class
class debug(class_registry.manage_state):
    """
    Debug parameters class for ISSM.

    This class encapsulates debugging and profiling parameters for the ISSM (Ice Sheet System Model) framework.
    It allows users to enable various debugging tools and profiling capabilities to analyze performance
    and troubleshoot issues in ice sheet model simulations.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    valgrind : int, default=0
        Use valgrind to debug (0 or 1).
    gprof : int, default=0
        Use GNU profiler to find out where the time is spent.
    profiling : int, default=0
        Enables profiling (memory, flops, time).

    Methods
    -------
    __init__(self, other=None)
        Initializes the debug parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the debug parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.debug = pyissm.param.debug()
    md.debug.valgrind = 1
    md.debug.profiling = 1
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.valgrind = 0
        self.gprof = 0
        self.profiling = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   debug parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'valgrind', 'use valgrind to debug (0 or 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'gprof', 'use gnu - profiler to find out where the time is spent'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'profiling', 'enables profiling (memory, flops, time)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - debug Class'
        return s

