import numpy as np
from . import param_utils
from . import class_registry
from .. import execute
from .. import utils

## ------------------------------------------------------
## inversion.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    """
    Default inversion parameters class for ISSM.

    This class defines the default parameters for the ISSM inversion process.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    iscontrol : int, default=0
        Is inversion activated? (0: no, 1: yes)
    incomplete_adjoint : int, default=1
        1: linear viscosity, 0: non-linear viscosity
    control_parameters : str, default='FrictionCoefficient'
        Control parameter(s) for inversion (e.g., 'FrictionCoefficient', 'MaterialsRheologyBbar')
    nsteps : int, default=20
        Number of optimization searches
    maxiter_per_step : ndarray, default=20*np.ones(nsteps)
        Maximum iterations during each optimization step
    cost_functions : str, default='List of cost functions'
        Type of response for each optimization step
    cost_functions_coefficients : ndarray, default=np.nan
        Coefficients applied to the misfit of each vertex and for each control parameter
    gradient_scaling : ndarray, default=50*np.ones((nsteps, 1))
        Scaling factor on gradient direction during optimization, for each optimization step
    cost_function_threshold : float, default=np.nan
        Misfit convergence criterion. Default is 1%, NaN if not applied
    min_parameters : float, default=np.nan
        Absolute minimum acceptable value of the inversed parameter on each vertex
    max_parameters : float, default=np.nan
        Absolute maximum acceptable value of the inversed parameter on each vertex
    step_threshold : ndarray, default=0.7*np.ones(nsteps)
        Decrease threshold for misfit, default is 30%
    vx_obs : ndarray, default=np.nan
        Observed velocity x component [m/yr]
    vy_obs : ndarray, default=np.nan
        Observed velocity y component [m/yr]
    vz_obs : ndarray, default=np.nan
        Observed velocity z component [m/yr]
    vel_obs : ndarray, default=np.nan
        Observed velocity magnitude [m/yr]
    thickness_obs : ndarray, default=np.nan
        Observed thickness [m]
    surface_obs : ndarray, default=np.nan
        Observed surface elevation [m]

    Methods
    -------
    __init__(self, other=None)
        Initializes the default inversion parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the default inversion parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.inversion = pyissm.param.inversion.default()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.iscontrol = 0
        self.incomplete_adjoint = 1
        self.control_parameters = 'FrictionCoefficient'
        self.nsteps = 20
        self.maxiter_per_step = 20 * np.ones(self.nsteps)
        self.cost_functions = 'List of cost functions' # Default = [101, ]
        self.cost_functions_coefficients = np.nan
        self.gradient_scaling = 50 * np.ones((self.nsteps, 1))
        self.cost_function_threshold = np.nan
        self.min_parameters = np.nan
        self.max_parameters = np.nan
        self.step_threshold = 0.7 * np.ones(self.nsteps)
        self.vx_obs = np.nan
        self.vy_obs = np.nan
        self.vz_obs = np.nan
        self.vel_obs = np.nan
        self.thickness_obs = np.nan
        self.surface_obs = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   inversion parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'iscontrol', 'is inversion activated?'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'incomplete_adjoint', '1: linear viscosity, 0: non - linear viscosity'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'control_parameters', 'ex: {''FrictionCoefficient''}, or {''MaterialsRheologyBbar''}'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'nsteps', 'number of optimization searches'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'cost_functions', 'indicate the type of response for each optimization step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'cost_functions_coefficients', 'cost_functions_coefficients applied to the misfit of each vertex and for each control_parameter'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'cost_function_threshold', 'misfit convergence criterion. Default is 1%, NaN if not applied'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'maxiter_per_step', 'maximum iterations during each optimization step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'gradient_scaling', 'scaling factor on gradient direction during optimization, for each optimization step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'step_threshold', 'decrease threshold for misfit, default is 30%'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'min_parameters', 'absolute minimum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'max_parameters', 'absolute maximum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vx_obs', 'observed velocity x component [m/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vy_obs', 'observed velocity y component [m/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vel_obs', 'observed velocity magnitude [m/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'thickness_obs', 'observed thickness [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'surface_obs', 'observed surface elevation [m]'))
        s += '{}\n'.format('Available cost functions:')
        s += '{}\n'.format('   101: SurfaceAbsVelMisfit')
        s += '{}\n'.format('   102: SurfaceRelVelMisfit')
        s += '{}\n'.format('   103: SurfaceLogVelMisfit')
        s += '{}\n'.format('   104: SurfaceLogVxVyMisfit')
        s += '{}\n'.format('   105: SurfaceAverageVelMisfit')
        s += '{}\n'.format('   201: ThicknessAbsMisfit')
        s += '{}\n'.format('   501: DragCoefficientAbsGradient')
        s += '{}\n'.format('   502: RheologyBbarAbsGradient')
        s += '{}\n'.format('   503: ThicknessAbsGradient')
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - inversion.default Class'
        return s
    
    # Marshall method for saving the inversion.default parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [inversion.default] parameters to a binary file.

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

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.inversion.type', data = 0, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'iscontrol', format = 'Boolean')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'incomplete_adjoint', format = 'Boolean')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'vel_obs', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts)

        ## Write conditional fields
        if self.iscontrol:
            
            ## Write DoubleMat fields
            execute.WriteData(fid, prefix, obj = self, fieldname = 'cost_functions_coefficients', format = 'DoubleMat', mattype = 1)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'gradient_scaling', format = 'DoubleMat', mattype = 3)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'min_parameters', format = 'DoubleMat', mattype = 3)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'max_parameters', format = 'DoubleMat', mattype = 3)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'step_threshold', format = 'DoubleMat', mattype = 3)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'vx_obs', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'vy_obs', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'vz_obs', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'thickness_obs', format = 'DoubleMat', mattype = 1)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'surface_obs', format = 'DoubleMat', mattype = 1)

            ## Write other fields
            execute.WriteData(fid, prefix, obj = self, fieldname = 'cost_function_threshold', format = 'Double')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'nsteps', format = 'Integer')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'maxiter_per_step', format = 'IntMat', mattype = 3)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'control_parameters', format = 'StringArray')
            execute.WriteData(fid, prefix, name = 'md.inversion.num_control_parameters', data = len(self.control_parameters), format = 'Integer')
            execute.WriteData(fid, prefix, name = 'md.inversion.cost_functions', data = param_utils.marshall_inversion_cost_functions(self.cost_functions), format = 'StringArray')
            execute.WriteData(fid, prefix, name = 'md.inversion.num_cost_functions', data = np.size(self.cost_functions), format = 'Integer')


## ------------------------------------------------------
## inversion.m1qn3
## ------------------------------------------------------
@class_registry.register_class
class m1qn3(class_registry.manage_state):
    """
    m1qn3 inversion parameters class for ISSM.

    This class defines the default parameters for the ISSM inversion process using the m1qn3 optimization algorithm.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    iscontrol : int, default=0
        Is inversion activated? (0: no, 1: yes)
    incomplete_adjoint : int, default=1
        1: linear viscosity, 0: non-linear viscosity
    control_parameters : str, default='FrictionCoefficient'
        Control parameter(s) for inversion (e.g., 'FrictionCoefficient', 'MaterialsRheologyBbar')
    control_scaling_factors : float or ndarray, default=1
        Order of magnitude of each control (useful for multi-parameter optimization)
    maxsteps : int, default=20
        Maximum number of iterations (gradient computation)
    maxiter : int, default=40
        Maximum number of function evaluations (forward run)
    dxmin : float, default=0.1
        Convergence criterion: two points less than dxmin from each other (sup-norm) are considered identical
    dfmin_frac : float, default=1.0
        Expected reduction of cost function during the first step (e.g., 0.3 = 30% reduction)
    gttol : float, default=1e4
        ||g(X)||/||g(X0)|| (g(X0): gradient at initial guess X0)
    cost_functions : int or list, default=101
        Type of response for each optimization step
    cost_functions_coefficients : ndarray, default=np.nan
        Coefficients applied to the misfit of each vertex and for each control parameter
    min_parameters : float, default=np.nan
        Absolute minimum acceptable value of the inversed parameter on each vertex
    max_parameters : float, default=np.nan
        Absolute maximum acceptable value of the inversed parameter on each vertex
    vx_obs : ndarray, default=np.nan
        Observed velocity x component [m/yr]
    vy_obs : ndarray, default=np.nan
        Observed velocity y component [m/yr]
    vz_obs : ndarray, default=np.nan
        Observed velocity z component [m/yr]
    vel_obs : ndarray, default=np.nan
        Observed velocity magnitude [m/yr]
    thickness_obs : ndarray, default=np.nan
        Observed thickness [m]
    surface_obs : ndarray, default=np.nan
        Observed surface elevation [m]

    Methods
    -------
    __init__(self, other=None)
        Initializes the m1qn3 inversion parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the m1qn3 inversion parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.inversion = pyissm.param.inversion.m1qn3()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.iscontrol = 0
        self.incomplete_adjoint = 1
        self.control_parameters = 'FrictionCoefficient'
        self.control_scaling_factors = 1
        self.maxsteps = 20
        self.maxiter = 40
        self.dxmin = 0.1
        self.dfmin_frac = 1.
        self.gttol = 1e4
        self.cost_functions = 101
        self.cost_functions_coefficients = np.nan
        self.min_parameters = np.nan
        self.max_parameters = np.nan
        self.vx_obs = np.nan
        self.vy_obs = np.nan
        self.vz_obs = np.nan
        self.vel_obs = np.nan
        self.thickness_obs = np.nan
        self.surface_obs = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   m1qn3 inversion parameters:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'iscontrol', 'is inversion activated?'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'incomplete_adjoint', '1: linear viscosity, 0: non - linear viscosity'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'control_parameters', 'ex: [''FrictionCoefficient''], or [''MaterialsRheologyBbar'']'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'control_scaling_factors', 'order of magnitude of each control (useful for multi - parameter optimization)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'maxsteps', 'maximum number of iterations (gradient computation)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'maxiter', 'maximum number of Function evaluation (forward run)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'dxmin', 'convergence criterion: two points less than dxmin from eachother (sup - norm) are considered identical'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'dfmin_frac', 'expected reduction of J during the first step (e.g., 0.3=30% reduction in cost function)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'gttol', '||g(X)||/||g(X0)|| (g(X0): gradient at initial guess X0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'cost_functions', 'indicate the type of response for each optimization step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'cost_functions_coefficients', 'cost_functions_coefficients applied to the misfit of each vertex and for each control_parameter'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'min_parameters', 'absolute minimum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'max_parameters', 'absolute maximum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vx_obs', 'observed velocity x component [m / yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vy_obs', 'observed velocity y component [m / yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vel_obs', 'observed velocity magnitude [m / yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'thickness_obs', 'observed thickness [m]'))
        s += '{}\n'.format('Available cost functions:')
        s += '{}\n'.format('   101: SurfaceAbsVelMisfit')
        s += '{}\n'.format('   102: SurfaceRelVelMisfit')
        s += '{}\n'.format('   103: SurfaceLogVelMisfit')
        s += '{}\n'.format('   104: SurfaceLogVxVyMisfit')
        s += '{}\n'.format('   105: SurfaceAverageVelMisfit')
        s += '{}\n'.format('   201: ThicknessAbsMisfit')
        s += '{}\n'.format('   501: DragCoefficientAbsGradient')
        s += '{}\n'.format('   502: RheologyBbarAbsGradient')
        s += '{}\n'.format('   503: ThicknessAbsGradient')
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - inversion.m1qn3 Class'
        return s

    # Marshall method for saving the inversion.m1qn3 parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [inversion.m1qn3] parameters to a binary file.

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

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.inversion.type', data = 2, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'iscontrol', format = 'Boolean')

        ## Write conditional fields
        if self.iscontrol:

            ## Write DoubleMat fields
            execute.WriteData(fid, prefix, obj = self, fieldname = 'cost_functions_coefficients', format = 'DoubleMat', mattype = 1)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'min_parameters', format = 'DoubleMat', mattype = 3)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'max_parameters', format = 'DoubleMat', mattype = 3)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'vx_obs', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'vy_obs', format = 'DoubleMat', mattype = 1, scale =  1. / md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'vz_obs', format = 'DoubleMat', mattype = 1, scale =  1. / md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'thickness_obs', format = 'DoubleMat', mattype = 1)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'surface_obs', format = 'DoubleMat', mattype = 1)

            ## Write other fields
            execute.WriteData(fid, prefix, obj = self, fieldname = 'incomplete_adjoint', format = 'Boolean')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'control_scaling_factors', format = 'DoubleMat', mattype = 3)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'maxsteps', format = 'Integer')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'maxiter', format = 'Integer')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'dxmin', format = 'Double')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'dfmin_frac', format = 'Double')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'gttol', format = 'Double')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'control_parameters', format = 'StringArray')
            execute.WriteData(fid, prefix, name = 'md.inversion.num_control_parameters', data = len(self.control_parameters), format = 'Integer')
            execute.WriteData(fid, prefix, name = 'md.inversion.cost_functions', data = param_utils.marshall_inversion_cost_functions(self.cost_functions), format = 'StringArray')
            execute.WriteData(fid, prefix, name = 'md.inversion.num_cost_functions', data = np.size(self.cost_functions), format = 'Integer')

## ------------------------------------------------------
## inversion.tao
## ------------------------------------------------------
@class_registry.register_class
class tao(class_registry.manage_state):
    """
    tao inversion parameters class for ISSM.

    This class defines the default parameters for the ISSM inversion process using the TAO optimization algorithms.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    iscontrol : int, default=0
        Is inversion activated? (0: no, 1: yes)
    incomplete_adjoint : int, default=1
        1: linear viscosity, 0: non-linear viscosity
    control_parameters : str, default='FrictionCoefficient'
        Control parameter(s) for inversion (e.g., 'FrictionCoefficient', 'MaterialsRheologyBbar')
    maxsteps : int, default=20
        Maximum number of iterations (gradient computation)
    maxiter : int, default=30
        Maximum number of function evaluations (forward run)
    fatol : float, default=0
        Absolute tolerance for cost function convergence (f(X) - f(X*))
    frtol : float, default=0
        Relative tolerance for cost function convergence (|f(X) - f(X*)| / |f(X*)|)
    gatol : float, default=0
        Absolute tolerance for gradient norm convergence (||g(X)||)
    grtol : float, default=0
        Relative tolerance for gradient norm convergence (||g(X)|| / |f(X)|)
    gttol : float, default=1e-4
        Tolerance for gradient norm relative to initial (||g(X)|| / ||g(X0)||)
    algorithm : str, default='blmvm'
        Minimization algorithm: 'blmvm', 'cg', 'lmvm'
    cost_functions : int or list, default=101
        Type of response for each optimization step
    cost_functions_coefficients : ndarray, default=np.nan
        Coefficients applied to the misfit of each vertex and for each control parameter
    min_parameters : float, default=np.nan
        Absolute minimum acceptable value of the inversed parameter on each vertex
    max_parameters : float, default=np.nan
        Absolute maximum acceptable value of the inversed parameter on each vertex
    vx_obs : ndarray, default=np.nan
        Observed velocity x component [m/yr]
    vy_obs : ndarray, default=np.nan
        Observed velocity y component [m/yr]
    vz_obs : ndarray, default=np.nan
        Observed velocity z component [m/yr]
    vel_obs : ndarray, default=np.nan
        Observed velocity magnitude [m/yr]
    thickness_obs : ndarray, default=np.nan
        Observed thickness [m]
    surface_obs : ndarray, default=np.nan
        Observed surface elevation [m]

    Methods
    -------
    __init__(self, other=None)
        Initializes the tao inversion parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the tao inversion parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.inversion = pyissm.param.inversion.tao()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.iscontrol = 0
        self.incomplete_adjoint = 1
        self.control_parameters = 'FrictionCoefficient'
        self.maxsteps = 20
        self.maxiter = 30
        self.fatol = 0
        self.frtol = 0
        self.gatol = 0
        self.grtol = 0
        self.gttol = 1e-4
        self.algorithm = self._get_algorithm()
        self.cost_functions = 101
        self.cost_functions_coefficients = np.nan
        self.min_parameters = np.nan
        self.max_parameters = np.nan
        self.vx_obs = np.nan
        self.vy_obs = np.nan
        self.vz_obs = np.nan
        self.vel_obs = np.nan
        self.thickness_obs = np.nan
        self.surface_obs = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   taoinversion parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'iscontrol', 'is inversion activated?'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'incomplete_adjoint', '1: linear viscosity, 0: non - linear viscosity'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'control_parameters', 'ex: {''FrictionCoefficient''}, or {''MaterialsRheologyBbar''}'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'maxsteps', 'maximum number of iterations (gradient computation)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'maxiter', 'maximum number of Function evaluation (forward run)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'fatol', 'convergence criterion: f(X) - f(X * ) (X: current iteration, X * : "true" solution, f: cost function)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'frtol', 'convergence criterion: |f(X) - f(X * )| / |f(X * )|'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'gatol', 'convergence criterion: ||g(X)|| (g: gradient of the cost function)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'grtol', 'convergence criterion: ||g(X)|| / |f(X)|'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'gttol', 'convergence criterion: ||g(X)|| / ||g(X0)|| (g(X0): gradient at initial guess X0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'algorithm', 'minimization algorithm: ''tao_blmvm'', ''tao_cg'', ''tao_lmvm'''))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'cost_functions', 'indicate the type of response for each optimization step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'cost_functions_coefficients', 'cost_functions_coefficients applied to the misfit of each vertex and for each control_parameter'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'min_parameters', 'absolute minimum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'max_parameters', 'absolute maximum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vx_obs', 'observed velocity x component [m / yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vy_obs', 'observed velocity y component [m / yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vel_obs', 'observed velocity magnitude [m / yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'thickness_obs', 'observed thickness [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'surface_obs', 'observed surface elevation [m]'))
        s += '{}\n'.format('Available cost functions:')
        s += '{}\n'.format('   101: SurfaceAbsVelMisfit')
        s += '{}\n'.format('   102: SurfaceRelVelMisfit')
        s += '{}\n'.format('   103: SurfaceLogVelMisfit')
        s += '{}\n'.format('   104: SurfaceLogVxVyMisfit')
        s += '{}\n'.format('   105: SurfaceAverageVelMisfit')
        s += '{}\n'.format('   201: ThicknessAbsMisfit')
        s += '{}\n'.format('   501: DragCoefficientAbsGradient')
        s += '{}\n'.format('   502: RheologyBbarAbsGradient')
        s += '{}\n'.format('   503: ThicknessAbsGradient')
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - inversion.tao Class'
        return s
    
    # Determine appropriate algorithm based on PETSc version
    def _get_algorithm(self):
        petsc_major = utils.wrappers.IssmConfig('_PETSC_MAJOR_')[0]
        petsc_minor = utils.wrappers.IssmConfig('_PETSC_MINOR_')[0]

        if petsc_major > 3 or (petsc_major == 3 and petsc_minor >= 5):
            algorithm = 'blmvm'
        else:
            algorithm = 'tao_blmvm'
        return algorithm
    
    # Marshall method for saving the inversion.tao parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [inversion.tao] parameters to a binary file.

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

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.inversion.type', data = 1, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'iscontrol', format = 'Boolean')

        ## Write conditional fields
        if self.iscontrol:

            ## Write Double fields
            fieldnames = ['fatol', 'frtol', 'gatol', 'grtol', 'gttol']
            for field in fieldnames:
                execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Double')

            ## Write DoubleMat fields
            execute.WriteData(fid, prefix, obj = self, fieldname = 'cost_functions_coefficients', format = 'DoubleMat', mattype = 1)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'min_parameters', format = 'DoubleMat', mattype = 3)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'max_parameters', format = 'DoubleMat', mattype = 3)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'vx_obs', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'vy_obs', format = 'DoubleMat', mattype = 1, scale =  1. / md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'vz_obs', format = 'DoubleMat', mattype = 1, scale =  1. / md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'thickness_obs', format = 'DoubleMat', mattype = 1)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'surface_obs', format = 'DoubleMat', mattype = 1)

            ## Write other fields
            execute.WriteData(fid, prefix, obj = self, fieldname = 'incomplete_adjoint', format = 'Boolean')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'maxsteps', format = 'Integer')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'maxiter', format = 'Integer')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'algorithm', format = 'String')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'control_parameters', format = 'StringArray')
            execute.WriteData(fid, prefix, name = 'md.inversion.num_control_parameters', data = len(self.control_parameters), format = 'Integer')
            execute.WriteData(fid, prefix, name = 'md.inversion.cost_functions', data = param_utils.marshall_inversion_cost_functions(self.cost_functions), format = 'StringArray')
            execute.WriteData(fid, prefix, name = 'md.inversion.num_cost_functions', data = np.size(self.cost_functions), format = 'Integer')
