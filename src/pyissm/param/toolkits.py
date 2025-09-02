import numpy as np
from . import class_registry
from . import param_utils
from .. import utils

@class_registry.register_class
class toolkits(class_registry.manage_state):
    """
    Toolkits parameters class for ISSM.

    This class encapsulates toolkit configuration parameters for the ISSM (Ice Sheet System Model) framework.
    Toolkits define solver options and configurations for different types of analyses, including linear solvers,
    preconditioners, and other computational tools used in ISSM simulations.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    DefaultToolkit : dict
        Default toolkit configuration (automatically set based on available libraries like MUMPS, GSL, PETSc).
    RecoveryAnalysis : dict
        Toolkit configuration for recovery mode analysis (same as DefaultToolkit by default).

    Methods
    -------
    __init__(self, other=None)
        Initializes the toolkits parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the toolkits parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.toolkits = pyissm.param.toolkits()
    # Toolkit configurations are automatically set based on available libraries
    # Additional analysis-specific toolkits can be added as needed
    """

    # Initialise with default parameters
    def __init__(self, other = None):

        ## Default toolkits
            ## Check for PETSc
        if utils.wrappers.IssmConfig('_HAVE_PETSC_')[0]:
            
            ## Check for toolkit
            if utils.wrappers.IssmConfig('_HAVE_MUMPS_')[0]:
                self.DefaultAnalysis = utils.config.mumps_options()
            else:
                self.DefaultAnalysis = utils.config.iluasm_options()

        else:
            if utils.wrappers.IssmConfig('_HAVE_MUMPS_')[0]:
                self.DefaultAnalysis = utils.config.issm_mumps_solver()
            elif utils.wrappers.IssmConfig('_HAVE_GSL_')[0]:
                self.DefaultAnalysis = utils.config.issm_gsl_solver()
            else:
                raise IOError(f'toolkits: need at least MUMPS or GSL to define ISSM solver type, no default solver assigned')

        ## Recovery mode (same as DefaultAnalysis by default)
        self.RecoveryAnalysis = self.DefaultAnalysis
                
        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = "List of toolkits options per analysis:\n\n"
        for analysis in list(vars(self).keys()):
            s += "{}\n".format(param_utils.fielddisplay(self, analysis, ''))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - toolkits Class'
        return s