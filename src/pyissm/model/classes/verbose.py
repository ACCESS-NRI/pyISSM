from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm.model import execute

@class_registry.register_class
class verbose(class_registry.manage_state):
    """
    Verbose parameters class for ISSM.

    This class encapsulates verbose parameters that control the level of output
    and logging in the ISSM (Ice Sheet System Model) framework. It allows users
    to enable or disable verbosity for different components of the model.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values
        in ``other`` differ from default values, they will override the default values.

    Attributes
    ----------
    mprocessor : :class:`int`, default=0
        Processor verbosity (0 or 1).
    module : :class:`int`, default=0
        Module verbosity (0 or 1).
    solution : :class:`int`, default=1
        Solution verbosity (0 or 1).
    solver : :class:`int`, default=0
        Solver verbosity (0 or 1).
    convergence : :class:`int`, default=0
        Convergence verbosity (0 or 1).
    control : :class:`int`, default=1
        Control verbosity (0 or 1).
    qmu : :class:`int`, default=1
        QMU (Dakota) verbosity (0 or 1).
    autodiff : :class:`int`, default=0
        Automatic differentiation verbosity (0 or 1).
    smb : :class:`int`, default=0
        Surface mass balance verbosity (0 or 1).

    Examples
    --------
    .. code-block:: python

        >>> md.verbose = pyissm.model.classes.verbose()
        >>> md.verbose.solution = 1
        >>> md.verbose.convergence = 1
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.mprocessor = 0
        self.module = 0
        self.solution = 1
        self.solver = 0
        self.convergence = 0
        self.control = 1
        self.qmu = 1
        self.autodiff = 0
        self.smb = 0

        # Inherit matching fields from provided class
        super().__init__(other)
    
        # Define repr
    def __repr__(self):
        s = '   verbose parameters:\n'

        s += '{}\n'.format(class_utils._field_display(self, 'mprocessor', 'processor verbosity'))
        s += '{}\n'.format(class_utils._field_display(self, 'module', 'module verbosity'))
        s += '{}\n'.format(class_utils._field_display(self, 'solution', 'solution verbosity'))
        s += '{}\n'.format(class_utils._field_display(self, 'solver', 'solver verbosity'))
        s += '{}\n'.format(class_utils._field_display(self, 'convergence', 'convergence verbosity'))
        s += '{}\n'.format(class_utils._field_display(self, 'control', 'control verbosity'))
        s += '{}\n'.format(class_utils._field_display(self, 'qmu', 'QMU verbosity'))
        s += '{}\n'.format(class_utils._field_display(self, 'autodiff', 'autodiff verbosity'))
        s += '{}\n'.format(class_utils._field_display(self, 'smb', 'SMB verbosity'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - verbose Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [verbose.verbose] parameters.

        Parameters
        ----------
        md : :class:`pyissm.model.Model`
            The model object to check.
        solution : :class:`str`
            The solution name to check.
        analyses : list of :class:`str`
            List of analyses to check consistency for.

        Returns
        -------
        md : :class:`pyissm.model.Model`
            The model object with any consistency errors noted.
        """
        # No checks
        return md
    
    def VerboseToBinary(self):
        """
        Convert current verbosity settings to integer bitmask.

        This method converts the boolean verbosity flags into a single integer
        where each bit represents a different verbosity option according to
        the field mapping defined in the class.

        Returns
        -------
        int
            Integer bitmask representing the current verbosity settings.
            Each bit corresponds to a verbosity flag as defined in _fields.

        Examples
        --------
        >>> verbose_obj = verbose()
        >>> verbose_obj.solution = True
        >>> verbose_obj.control = True
        >>> binary_value = verbose_obj.VerboseToBinary()
        >>> print(binary_value)  # Will print 36 (4 + 32)
        """
        
        binary = 0
        if self.mprocessor:
            binary = binary | 1
        if self.module:
            binary = binary | 2
        if self.solution:
            binary = binary | 4
        if self.solver:
            binary = binary | 8
        if self.convergence:
            binary = binary | 16
        if self.control:
            binary = binary | 32
        if self.qmu:
            binary = binary | 64
        if self.autodiff:
            binary = binary | 128
        if self.smb:
            binary = binary | 256

        return binary

    def deactivate_all(self):
        """
        Deactivate all verbose model components.
        """

        self.mprocessor = 0
        self.module = 0
        self.solution = 0
        self.solver = 0
        self.convergence = 0
        self.control = 0
        self.qmu = 0
        self.autodiff = 0
        self.smb = 0

        return self

    # Marshall method for saving the verbose parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [verbose.verbose] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
            ISSM model object needed in some cases.

        Returns
        -------
        None
        """

        ## Write fields
        execute._write_model_field(fid, prefix, name = 'md.verbose', data = self.VerboseToBinary(), format = 'Integer')