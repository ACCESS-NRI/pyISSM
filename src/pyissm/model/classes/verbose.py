from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm.model import execute

@class_registry.register_class
class verbose(class_registry.manage_state):
    """
    Verbose parameters class for ISSM.

    This class encapsulates verbose parameters that control the level of output and logging in the ISSM (Ice Sheet System Model) framework.
    It stores settings for logging frequency, output file paths, and other related options.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    surface : ndarray, default=nan
        Ice upper surface elevation [m].
    thickness : ndarray, default=nan
        Ice thickness [m].
    base : ndarray, default=nan
        Ice base elevation [m].
    bed : ndarray, default=nan
        Bed elevation [m].
    hydrostatic_ratio : float, default=nan
        Hydrostatic ratio for floating ice.

    Methods
    -------
    __init__(self, other=None)
        Initializes the geometry parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the geometry parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.geometry = pyissm.model.classes.geometry()
    md.geometry.surface = surface_elevation
    md.geometry.thickness = ice_thickness
    md.geometry.bed = bed_elevation
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
        Marshall [verbose] parameters to a binary file.

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
        execute.WriteData(fid, prefix, name = 'md.verbose', data = self.VerboseToBinary(), format = 'Integer')