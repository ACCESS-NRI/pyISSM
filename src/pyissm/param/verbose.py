import numpy as np
from . import param_utils
from . import class_registry
from .. import execute

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

    Examples
    --------
    md.geometry = pyissm.param.geometry()
    md.geometry.surface = surface_elevation
    md.geometry.thickness = ice_thickness
    md.geometry.bed = bed_elevation
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.mprocessor = False
        self.module = False
        self.solution = False
        self.solver = False
        self.convergence = False
        self.control = False
        self.qmu = False
        self.autodiff = False
        self.smb = False

        # Inherit matching fields from provided class
        super().__init__(other)
    
        # Define repr
    def __repr__(self):
        s = '   verbose parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'mprocessor', 'processor verbosity'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'module', 'module verbosity'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'solution', 'solution verbosity'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'solver', 'solver verbosity'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'convergence', 'convergence verbosity'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'control', 'control verbosity'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'qmu', 'QMU verbosity'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'autodiff', 'autodiff verbosity'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'smb', 'SMB verbosity'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - verbose Class'
        return s
    
    def VerboseToBinary(self):  # {{{
        #BEGINVERB2BIN
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
        #ENDVERB2BIN

        return binary

    # Marshall method for saving the verbose parameters
    def marshall_class(self, prefix, md, fid):
        """
        Marshall the verbose parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.

        Returns
        -------
        None
        """

        execute.WriteData(fid, prefix, name = 'md.verbose', data = self.VerboseToBinary(), format = 'Integer')