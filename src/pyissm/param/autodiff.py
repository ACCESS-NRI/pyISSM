import numpy as np
from . import param_utils
from . import class_registry
from .. import execute

@class_registry.register_class
class autodiff(class_registry.manage_state):
    """
    Automatic differentiation parameters class for ISSM.

    This class encapsulates parameters for automatic differentiation (AD) functionality in the ISSM (Ice Sheet System Model) framework.
    It allows users to configure AD settings including dependent and independent variables, memory buffer sizes, 
    and optimization parameters for sensitivity analysis and gradient-based optimization.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    isautodiff : float, default=0.0
        Indicates if automatic differentiation is activated.
    dependents : str, default='List dependents'
        List of dependent variables for AD.
    independents : str, default='List independents'
        List of independent variables for AD.
    driver : str, default='fos_forward'
        ADOLC driver ('fos_forward' or 'fov_forward').
    obufsize : float, default=524288
        Number of operations per buffer (== OBUFSIZE in usrparms.h).
    lbufsize : float, default=524288
        Number of locations per buffer (== LBUFSIZE in usrparms.h).
    cbufsize : float, default=524288
        Number of values per buffer (== CBUFSIZE in usrparms.h).
    tbufsize : float, default=524288
        Number of taylors per buffer (<=TBUFSIZE in usrparms.h).
    gcTriggerMaxSize : float, default=65536
        Free location block sorting/consolidation triggered if allocated locations exceed this value.
    gcTriggerRatio : float, default=2.0
        Free location block sorting/consolidation triggered if the ratio between allocated and used locations exceeds this value.
    tapeAlloc : float, default=15000000
        Iteration count of a priori memory allocation of the AD tape.
    outputTapeMemory : float, default=0.0
        Write AD tape memory statistics to file ad_mem.dat.
    outputTime : float, default=0.0
        Write AD recording and evaluation times to file ad_time.dat.
    enablePreaccumulation : float, default=0.0
        Enable CoDiPack preaccumulation in augmented places.

    Methods
    -------
    __init__(self, other=None)
        Initializes the autodiff parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the autodiff parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.autodiff = pyissm.param.autodiff()
    md.autodiff.isautodiff = 1
    md.autodiff.dependents = ['Vel']
    md.autodiff.independents = ['MaterialsRheologyBbar']
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.isautodiff = 0.
        self.dependents = 'List dependents'
        self.independents = 'List independents'
        self.driver = 'fos_forward'
        self.obufsize = 524288
        self.lbufsize = 524288
        self.cbufsize = 524288
        self.tbufsize = 524288
        self.gcTriggerMaxSize = 65536
        self.gcTriggerRatio = 2.0
        self.tapeAlloc = 15000000
        self.outputTapeMemory = 0.
        self.outputTime = 0.
        self.enablePreaccumulation = 0.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):  # {{{
        s = '      automatic differentiation parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'isautodiff', "indicates if the automatic differentiation is activated"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'dependents', "list of dependent variables"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'independents', "list of independent variables"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'driver', "ADOLC driver ('fos_forward' or 'fov_forward')"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'obufsize', "Number of operations per buffer (== OBUFSIZE in usrparms.h)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'lbufsize', "Number of locations per buffer (== LBUFSIZE in usrparms.h)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'cbufsize', "Number of values per buffer (== CBUFSIZE in usrparms.h)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'tbufsize', "Number of taylors per buffer (<=TBUFSIZE in usrparms.h)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'gcTriggerRatio', "free location block sorting / consolidation triggered if the ratio between allocated and used locations exceeds gcTriggerRatio"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'gcTriggerMaxSize', "free location block sorting / consolidation triggered if the allocated locations exceed gcTriggerMaxSize)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'tapeAlloc', 'Iteration count of a priori memory allocation of the AD tape'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'outputTapeMemory', 'Write AD tape memory statistics to file ad_mem.dat'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'outputTime', 'Write AD recording and evaluation times to file ad_time.dat'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'enablePreaccumulation', 'Enable CoDiPack preaccumulation in augmented places'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - autodiff Class'
        return s

    # Marshall method for saving the autodiff parameters
    def marshall_class(self, prefix, md, fid):
        """
        Marshall the autodiff parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.

        Returns
        -------
        None
        """
        execute.WriteData(fid, prefix, obj = self, fieldname = 'isautodiff', format = 'Boolean')