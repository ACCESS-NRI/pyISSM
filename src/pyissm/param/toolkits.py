import warnings
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

        ## If wrappers are installed, set default toolkits based on available libraries
        if utils.wrappers.check_wrappers_installed():
            
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

        else:
            ## If wrappers are not installed, return empty toolkits class
            warnings.warn('pyissm.param.toolkits: Python wrappers not installed. Unable to automatically determine toolkit information.\n'
                          'Returning empty toolkits class to be defined manually.')
            self.DefaultAnalysis = {'toolkit': '',         
                                    'mat_type': '',
                                    'vec_type': '',
                                    'solver_type': ''
                                    }
            self.RecoveryAnalysis = self.DefaultAnalysis

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
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):

        supported_analyses = param_utils.supported_analyses()

        analyses = list(vars(self).keys())
        for analysis in analyses:
            if analysis not in supported_analyses:
                md.checkmessage('md.toolkits.{} not supported yet'.format(analysis))

            if not getattr(self, analysis):
                md.checkmessage('md.toolkits.{} is empty'.format(analysis))

        return md
    
    def write_toolkits_file(self, filename):
        """
        Build a PETSc compatible options file from the toolkits model field.

        This method creates a PETSc compatible options file from the current toolkits configuration.
        The generated file can be used with both 'petsc' and 'issm' toolkit types for solver
        configuration in ISSM simulations.

        Parameters
        ----------
        filename : str
            Path and name of the output file to write the toolkit options to.

        Raises
        ------
        IOError
            If the specified file cannot be opened for writing.
        TypeError
            If any option value is not a supported type (bool, int, float, or str).

        Examples
        --------
        >>> md.toolkits.write_toolkits_file('solver_options.toolkits')
        """

        # Open file for writing
        try:
            fid = open(filename, 'w')
        except IOError as e:
            raise IOError(f'write_toolkits_file: could not open {filename} for writing: {e}')

        # Write header
        fid.write('{}{}{}\n'.format('%Toolkits options file: ', filename, ' written from Python toolkits array'))

        # Start writing options
        for analysis in list(vars(self).keys()):
            options = getattr(self, analysis)

            # Write analysis
            ## NOTE: Append a + to recognize it's an analysis enum
            fid.write('\n+{}\n'.format(analysis))

            # Write options
            for optionname, optionvalue in list(options.items()):
                
                ## If the option has no value, write the option name only
                if not optionvalue:
                    fid.write('-{}\n'.format(optionname))
                ## If the option has a name and a value, write both
                else:
                    ## optionvalue can be string or scalar
                    if isinstance(optionvalue, (bool, int, float, str)):
                        fid.write('-{} {}\n'.format(optionname, optionvalue))
                    else:
                        raise TypeError(f'write_toolkits_file: option {optionname} is not formatted correctly. Values must be string or scalar.')

        fid.close()

        