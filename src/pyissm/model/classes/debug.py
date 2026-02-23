from pyissm.model.classes import class_utils, class_registry
from pyissm.model import execute

@class_registry.register_class
class debug(class_registry.manage_state):
    """
    Debug class for ISSM.

    This class contains debugging and profiling parameters for the ISSM framework.
    It allows users to enable various debugging tools and profiling capabilities
    to analyze performance and troubleshoot issues in ice sheet model simulations.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    valgrind : :class:`int`, default=0
        Use valgrind to debug (0 or 1).
    gprof : :class:`int`, default=0
        Use GNU profiler to find out where the time is spent.
    profiling : :class:`int`, default=0
        Enables profiling (memory, flops, time).

    Examples
    --------
    .. code-block:: python
        >>> md.debug = pyissm.model.classes.debug()
        >>> md.debug.valgrind = 1
        >>> md.debug.profiling = 1
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

        s += '{}\n'.format(class_utils.fielddisplay(self, 'valgrind', 'use valgrind to debug (0 or 1)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'gprof', 'use gnu - profiler to find out where the time is spent'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'profiling', 'enables profiling (memory, flops, time)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - debug Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [debug] parameters.

        Parameters
        ----------
        md : :class:`pyissm.model.Model`
            The model object to check.
        solution : :class:`pyissm.model.solution`
            The solution object to check.
        analyses : list of :class:`str`
            List of analyses to check consistency for.

        Returns 
        -------
        md : :class:`pyissm.model.Model`
            The model object with any consistency errors noted.
        """

        # No specific checks for debug parameters
        return md
    
    # Marshall method for saving the debug parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [debug] parameters to a binary file.

        Parameters
        ----------
        fid : :class:`file object`
            The file object to write the binary data to.
        prefix : :class:`str`
            Prefix string used for data identification in the binary file.
        md : :class:`pyissm.model.Model`, optional
            ISSM model object needed in some cases.
            
        Returns
        -------
        None
        """

        ## Write single field to file
        execute.WriteData(fid, prefix, obj = self, fieldname = 'profiling', format = 'Boolean')