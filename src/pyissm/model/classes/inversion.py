import numpy as np
import warnings
from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm.model import execute, mesh
from pyissm import tools

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
    md.inversion = pyissm.model.classes.inversion.default()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.iscontrol = 0
        self.incomplete_adjoint = 1
        self.control_parameters = ['FrictionCoefficient']
        self.nsteps = 20
        self.maxiter_per_step = 20 * np.ones(self.nsteps)
        self.cost_functions = [101]
        self.cost_functions_coefficients = np.nan
        self.gradient_scaling = 50 * np.ones(self.nsteps)
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

        s += '{}\n'.format(class_utils.fielddisplay(self, 'iscontrol', 'is inversion activated?'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'incomplete_adjoint', '1: linear viscosity, 0: non - linear viscosity'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'control_parameters', 'ex: {''FrictionCoefficient''}, or {''MaterialsRheologyBbar''}'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'nsteps', 'number of optimization searches'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'cost_functions', 'indicate the type of response for each optimization step'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'cost_functions_coefficients', 'cost_functions_coefficients applied to the misfit of each vertex and for each control_parameter'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'cost_function_threshold', 'misfit convergence criterion. Default is 1%, NaN if not applied'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'maxiter_per_step', 'maximum iterations during each optimization step'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'gradient_scaling', 'scaling factor on gradient direction during optimization, for each optimization step'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'step_threshold', 'decrease threshold for misfit, default is 30%'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'min_parameters', 'absolute minimum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'max_parameters', 'absolute maximum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'vx_obs', 'observed velocity x component [m/yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'vy_obs', 'observed velocity y component [m/yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'vel_obs', 'observed velocity magnitude [m/yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'thickness_obs', 'observed thickness [m]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'surface_obs', 'observed surface elevation [m]'))
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
    
    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude inversion.default fields to 3D
        """
        self.vx_obs = mesh.project_3d(md, vector = self.vx_obs, type = 'node')
        self.vy_obs = mesh.project_3d(md, vector = self.vy_obs, type = 'node')
        self.vel_obs = mesh.project_3d(md, vector = self.vel_obs, type = 'node')
        self.thickness_obs = mesh.project_3d(md, vector = self.thickness_obs, type = 'node')
        if not np.any(np.isnan(self.cost_functions_coefficients)):
            self.cost_functions_coefficients = mesh.project_3d(md, vector = self.cost_functions_coefficients, type = 'node')
        if not np.any(np.isnan(self.min_parameters)):
            self.min_parameters = mesh.project_3d(md, vector = self.min_parameters, type = 'node')
        if not np.any(np.isnan(self.max_parameters)):
            self.max_parameters = mesh.project_3d(md, vector = self.max_parameters, type = 'node')
            
        return self

    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        # Early return if inversion disabled
        if not self.iscontrol:
            return md

        num_controls = np.size(md.inversion.control_parameters)
        num_costfunc = np.size(md.inversion.cost_functions)

        class_utils.check_field(md, fieldname = 'inversion.iscontrol', values = [0, 1])
        class_utils.check_field(md, fieldname = 'inversion.incomplete_adjoint', values = [0, 1])
        class_utils.check_field(md, fieldname = 'inversion.control_parameters', cell = True, values = class_utils.supported_inversion_control_parameters())
        class_utils.check_field(md, fieldname = 'inversion.nsteps', scalar = True, ge = 0)
        class_utils.check_field(md, fieldname = 'inversion.maxiter_per_step', size = (md.inversion.nsteps, ), ge = 0)
        class_utils.check_field(md, fieldname = 'inversion.step_threshold', size = (md.inversion.nsteps, ))
        class_utils.check_field(md, fieldname = 'inversion.cost_functions', size = (num_costfunc, ), values = class_utils.supported_inversion_cost_functions())
        if num_costfunc == 1:
            md.inversion.cost_functions_coefficients = np.squeeze(md.inversion.cost_functions_coefficients)
            class_utils.check_field(md, fieldname = 'inversion.cost_functions_coefficients', size = (md.mesh.numberofvertices, ), ge = 0)
        else:
            class_utils.check_field(md, fieldname = 'inversion.cost_functions_coefficients', size = (md.mesh.numberofvertices, num_costfunc), ge = 0)

        if num_controls == 1:
            md.inversion.gradient_scaling = np.squeeze(md.inversion.gradient_scaling)
            md.inversion.min_parameters = np.squeeze(md.inversion.min_parameters)
            md.inversion.max_parameters = np.squeeze(md.inversion.max_parameters)
            class_utils.check_field(md, fieldname = 'inversion.gradient_scaling', size = (md.inversion.nsteps, ))
            class_utils.check_field(md, fieldname = 'inversion.min_parameters', size = (md.mesh.numberofvertices, ))
            class_utils.check_field(md, fieldname = 'inversion.max_parameters', size = (md.mesh.numberofvertices, ))
        else:
            class_utils.check_field(md, fieldname = 'inversion.gradient_scaling', size = (md.inversion.nsteps, num_controls))
            class_utils.check_field(md, fieldname = 'inversion.min_parameters', size = (md.mesh.numberofvertices, num_controls))
            class_utils.check_field(md, fieldname = 'inversion.max_parameters', size = (md.mesh.numberofvertices, num_controls))

        # Only SSA, MMOLHO, L1L2, HO and FS are supported right now
        if solution == 'StressbalanceSolution':
            if not (md.flowequation.isSSA or md.flowequation.isMOLHO or md.flowequation.isHO or md.flowequation.isFS or md.flowequation.isL1L2):
                md.check_message("inversion can only be performed for SSA, MOLHO, HO or FS ice flow models")
        
        # Balancethicknesssolution
        if solution == 'BalancethicknessSolution':
            class_utils.check_field(md, fieldname = 'inversion.thickness_obs', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        elif solution == 'BalancethicknessSoftSolution':
            class_utils.check_field(md, fieldname = 'inversion.thickness_obs', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        else:
            class_utils.check_field(md, fieldname = 'inversion.vx_obs', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = 'inversion.vy_obs', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        return md
    
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
            execute.WriteData(fid, prefix, name = 'md.inversion.cost_functions', data = class_utils.marshall_inversion_cost_functions(self.cost_functions), format = 'StringArray')
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
    md.inversion = pyissm.model.classes.inversion.m1qn3()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.iscontrol = 0
        self.incomplete_adjoint = 1
        self.control_parameters = ['FrictionCoefficient']
        self.control_scaling_factors = 1
        self.maxsteps = 20
        self.maxiter = 40
        self.dxmin = 0.1
        self.dfmin_frac = 1.
        self.gttol = 1e-4
        self.cost_functions = [101]
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
        s += '{}\n'.format(class_utils.fielddisplay(self, 'iscontrol', 'is inversion activated?'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'incomplete_adjoint', '1: linear viscosity, 0: non - linear viscosity'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'control_parameters', 'ex: [''FrictionCoefficient''], or [''MaterialsRheologyBbar'']'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'control_scaling_factors', 'order of magnitude of each control (useful for multi - parameter optimization)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'maxsteps', 'maximum number of iterations (gradient computation)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'maxiter', 'maximum number of Function evaluation (forward run)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'dxmin', 'convergence criterion: two points less than dxmin from eachother (sup - norm) are considered identical'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'dfmin_frac', 'expected reduction of J during the first step (e.g., 0.3=30% reduction in cost function)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'gttol', '||g(X)||/||g(X0)|| (g(X0): gradient at initial guess X0)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'cost_functions', 'indicate the type of response for each optimization step'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'cost_functions_coefficients', 'cost_functions_coefficients applied to the misfit of each vertex and for each control_parameter'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'min_parameters', 'absolute minimum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'max_parameters', 'absolute maximum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'vx_obs', 'observed velocity x component [m / yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'vy_obs', 'observed velocity y component [m / yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'vel_obs', 'observed velocity magnitude [m / yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'thickness_obs', 'observed thickness [m]'))
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
    
    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude inversion.m1qn3 fields to 3D
        """
        self.vx_obs = mesh.project_3d(md, vector = self.vx_obs, type = 'node')
        self.vy_obs = mesh.project_3d(md, vector = self.vy_obs, type = 'node')
        self.vel_obs = mesh.project_3d(md, vector = self.vel_obs, type = 'node')
        self.thickness_obs = mesh.project_3d(md, vector = self.thickness_obs, type = 'node')
        if not np.any(np.isnan(self.cost_functions_coefficients)):
            self.cost_functions_coefficients = mesh.project_3d(md, vector = self.cost_functions_coefficients, type = 'node')
        if not np.any(np.isnan(self.min_parameters)):
            self.min_parameters = mesh.project_3d(md, vector = self.min_parameters, type = 'node')
        if not np.any(np.isnan(self.max_parameters)):
            self.max_parameters = mesh.project_3d(md, vector = self.max_parameters, type = 'node')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        # Early return if inversion disabled
        if not self.iscontrol:
            return md

        num_controls = np.size(md.inversion.control_parameters)
        num_costfunc = np.size(md.inversion.cost_functions)

        class_utils.check_field(md, fieldname = 'inversion.iscontrol', values = [0, 1])
        class_utils.check_field(md, fieldname = 'inversion.incomplete_adjoint', values = [0, 1])
        class_utils.check_field(md, fieldname = 'inversion.control_parameters', cell = True, values = class_utils.supported_inversion_control_parameters())
        class_utils.check_field(md, fieldname = 'inversion.control_scaling_factors', size = (num_controls, ), gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'inversion.maxsteps', scalar = True, ge = 0)
        class_utils.check_field(md, fieldname = 'inversion.maxiter', scalar = True, ge = 0)
        class_utils.check_field(md, fieldname = 'inversion.dxmin', scalar = True, gt = 0.)
        class_utils.check_field(md, fieldname = 'inversion.dfmin_frac', scalar = True, ge = 0., le = 1.)
        class_utils.check_field(md, fieldname = 'inversion.gttol', scalar = True, gt = 0.)
        class_utils.check_field(md, fieldname = 'inversion.cost_functions', size = (num_costfunc, ), values = class_utils.supported_inversion_cost_functions())
        if num_costfunc == 1:
            md.inversion.cost_functions_coefficients = np.squeeze(md.inversion.cost_functions_coefficients)
            class_utils.check_field(md, fieldname = 'inversion.cost_functions_coefficients', size = (md.mesh.numberofvertices, ), ge = 0)
        else:
            class_utils.check_field(md, fieldname = 'inversion.cost_functions_coefficients', size = (md.mesh.numberofvertices, num_costfunc), ge = 0)

        if num_controls == 1:
            md.inversion.min_parameters = np.squeeze(md.inversion.min_parameters)
            md.inversion.max_parameters = np.squeeze(md.inversion.max_parameters)
            class_utils.check_field(md, fieldname = 'inversion.min_parameters', size = (md.mesh.numberofvertices, ))
            class_utils.check_field(md, fieldname = 'inversion.max_parameters', size = (md.mesh.numberofvertices, ))
        else:
            class_utils.check_field(md, fieldname = 'inversion.min_parameters', size = (md.mesh.numberofvertices, num_controls))
            class_utils.check_field(md, fieldname = 'inversion.max_parameters', size = (md.mesh.numberofvertices, num_controls))

        if solution == 'BalancethicknessSolution':
            class_utils.check_field(md, fieldname = 'inversion.thickness_obs', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = 'inversion.surface_obs', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        elif solution == 'BalancethicknessSoftSolution':
            class_utils.check_field(md, fieldname = 'inversion.thickness_obs', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        else:
            class_utils.check_field(md, fieldname = 'inversion.vx_obs', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = 'inversion.vy_obs', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        
        return md

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
            execute.WriteData(fid, prefix, obj = self, fieldname = 'vel_obs', format = 'DoubleMat', mattype = 1, scale =  1. / md.constants.yts)

            ## Write conditional mattype
            if np.size(self.thickness_obs) == md.mesh.numberofelements:
                mattype = 2
            else:
                mattype = 1
            execute.WriteData(fid, prefix, obj = self, fieldname = 'thickness_obs', format = 'DoubleMat', mattype = mattype)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'surface_obs', format = 'DoubleMat', mattype = mattype)

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
            execute.WriteData(fid, prefix, name = 'md.inversion.cost_functions', data = class_utils.marshall_inversion_cost_functions(self.cost_functions), format = 'StringArray')
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
    md.inversion = pyissm.model.classes.inversion.tao()
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

        s += '{}\n'.format(class_utils.fielddisplay(self, 'iscontrol', 'is inversion activated?'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'incomplete_adjoint', '1: linear viscosity, 0: non - linear viscosity'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'control_parameters', 'ex: {''FrictionCoefficient''}, or {''MaterialsRheologyBbar''}'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'maxsteps', 'maximum number of iterations (gradient computation)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'maxiter', 'maximum number of Function evaluation (forward run)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'fatol', 'convergence criterion: f(X) - f(X * ) (X: current iteration, X * : "true" solution, f: cost function)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'frtol', 'convergence criterion: |f(X) - f(X * )| / |f(X * )|'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'gatol', 'convergence criterion: ||g(X)|| (g: gradient of the cost function)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'grtol', 'convergence criterion: ||g(X)|| / |f(X)|'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'gttol', 'convergence criterion: ||g(X)|| / ||g(X0)|| (g(X0): gradient at initial guess X0)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'algorithm', 'minimization algorithm: ''tao_blmvm'', ''tao_cg'', ''tao_lmvm'''))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'cost_functions', 'indicate the type of response for each optimization step'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'cost_functions_coefficients', 'cost_functions_coefficients applied to the misfit of each vertex and for each control_parameter'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'min_parameters', 'absolute minimum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'max_parameters', 'absolute maximum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'vx_obs', 'observed velocity x component [m / yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'vy_obs', 'observed velocity y component [m / yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'vel_obs', 'observed velocity magnitude [m / yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'thickness_obs', 'observed thickness [m]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'surface_obs', 'observed surface elevation [m]'))
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

    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude inversion.tao fields to 3D
        """
        self.vx_obs = mesh.project_3d(md, vector = self.vx_obs, type = 'node')
        self.vy_obs = mesh.project_3d(md, vector = self.vy_obs, type = 'node')
        self.vel_obs = mesh.project_3d(md, vector = self.vel_obs, type = 'node')
        self.thickness_obs = mesh.project_3d(md, vector = self.thickness_obs, type = 'node')
        
        if np.size(self.cost_functions_coefficients) > 1:
            self.cost_functions_coefficients = mesh.project_3d(md, vector = self.cost_functions_coefficients, type = 'node')

        if np.size(self.min_parameters) > 1:
            self.min_parameters = mesh.project_3d(md, vector = self.min_parameters, type = 'node')

        if np.size(self.max_parameters) > 1:
            self.max_parameters = mesh.project_3d(md, vector = self.max_parameters, type = 'node')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        
        # Early return if inversion disabled
        if not self.iscontrol:
            return md
        
        if tools.wrappers.check_wrappers_installed():                
            if not tools.wrappers.IssmConfig('_HAVE_TAO_')[0]:
                md = md.check_message('TAO has not been installed, ISSM needs to be reconfigured and recompiled with TAO')
        else:
            warnings.warn('pyissm.model.classes.inversion.tao.check_consistency: Python wrappers not installed. Unable to check for TAO installation.\n'
                          'Proceeding without TAO installation check.')
            
        num_controls = np.size(md.inversion.control_parameters)
        num_costfunc = np.size(md.inversion.cost_functions)

        class_utils.check_field(md, fieldname = 'inversion.iscontrol', values = [0, 1])
        class_utils.check_field(md, fieldname = 'inversion.incomplete_adjoint', values = [0, 1])
        class_utils.check_field(md, fieldname = 'inversion.control_parameters', cell = True, values = class_utils.supported_inversion_control_parameters())
        class_utils.check_field(md, fieldname = 'inversion.maxsteps', scalar = True, ge = 0)
        class_utils.check_field(md, fieldname = 'inversion.maxiter', scalar = True, ge = 0)
        class_utils.check_field(md, fieldname = 'inversion.fatol', scalar = True, ge = 0)
        class_utils.check_field(md, fieldname = 'inversion.frtol', scalar = True, ge = 0)
        class_utils.check_field(md, fieldname = 'inversion.gatol', scalar = True, ge = 0)
        class_utils.check_field(md, fieldname = 'inversion.grtol', scalar = True, ge = 0)
        class_utils.check_field(md, fieldname = 'inversion.gttol', scalar = True, ge = 0)

        if tools.wrappers.check_wrappers_installed():
            PETSCMAJOR = tools.wrappers.IssmConfig('_PETSC_MAJOR_')[0]
            PETSCMINOR = tools.wrappers.IssmConfig('_PETSC_MINOR_')[0]
            if(PETSCMAJOR > 3 or (PETSCMAJOR == 3 and PETSCMINOR >= 5)):
                class_utils.check_field(md, fieldname = 'inversion.algorithm', values = ['blmvm', 'cg', 'lmvm'])
            else:
                class_utils.check_field(md, fieldname = 'inversion.algorithm', values = ['tao_blmvm', 'tao_cg', 'tao_lmvm'])
        else:
            warnings.warn('pyissm.model.classes.inversion.tao.check_consistency: Python wrappers not installed. Unable to check PETSc version for algorithm validation.\n'
                          'Proceeding without algorithm validation.')

        class_utils.check_field(md, fieldname = 'inversion.cost_functions', size = (num_costfunc, ), values = class_utils.supported_inversion_cost_functions())
        class_utils.check_field(md, fieldname = 'inversion.cost_functions_coefficients', size = (md.mesh.numberofvertices, num_costfunc), ge = 0)
        class_utils.check_field(md, fieldname = 'inversion.min_parameters', size = (md.mesh.numberofvertices, num_controls))
        class_utils.check_field(md, fieldname = 'inversion.max_parameters', size = (md.mesh.numberofvertices, num_controls))

        if solution == 'BalancethicknessSolution':
            class_utils.check_field(md, fieldname = 'inversion.thickness_obs', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        elif solution == 'BalancethicknessSoftSolution':
            class_utils.check_field(md, fieldname = 'inversion.thickness_obs', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        else:
            class_utils.check_field(md, fieldname = 'inversion.vx_obs', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = 'inversion.vy_obs', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)

        return md
    
    # Determine appropriate algorithm based on PETSc version
    def _get_algorithm(self):

        ## If python wrappers are installed, check PETSc version
        if tools.wrappers.check_wrappers_installed():
            petsc_major = tools.wrappers.IssmConfig('_PETSC_MAJOR_')[0]
            petsc_minor = tools.wrappers.IssmConfig('_PETSC_MINOR_')[0]

            if petsc_major > 3 or (petsc_major == 3 and petsc_minor >= 5):
                algorithm = 'blmvm'
            else:
                algorithm = 'tao_blmvm'
        else:
            ## If python wrappers are not installed, return empty string and let ISSM decide
            warnings.warn('pyissm.model.classes.inversion.tao: Python wrappers not installed. Unable to automatically determine algorithm.\n'
                          'Returning empty algorithm to be defined manually.')
            algorithm = ''
        
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
            execute.WriteData(fid, prefix, name = 'md.inversion.cost_functions', data = class_utils.marshall_inversion_cost_functions(self.cost_functions), format = 'StringArray')
            execute.WriteData(fid, prefix, name = 'md.inversion.num_cost_functions', data = np.size(self.cost_functions), format = 'Integer')

@class_registry.register_class
class adm1qn3(class_registry.manage_state):
    """
    ADM1QN3 inversion parameters for ISSM.

    This class configures the **ADM1QN3** inversion driver (ISSM inversion type = 4).
    It is the Python equivalent of the MATLAB ``adm1qn3inversion`` class, but follows
    the same pyISSM conventions as other inversion classes (e.g. ``inversion.m1qn3``,
    ``inversion.tao``):

    - construction via defaults + optional inheritance from ``other``
    - ``check_consistency`` validation
    - ``marshall_class`` serialization to ISSM binary format
    - a no-op ``extrude`` for interface consistency

    ADM1QN3 provides a quasi-Newton style optimization loop with tunable iteration
    limits and convergence thresholds. This class only stores the *optimizer-level*
    settings; model-dependent inversion fields (e.g. observations, cost functions,
    control parameters) are handled by other inversion types/classes.

    Parameters
    ----------
    other : object, optional
        Any object with matching attribute names. If provided, its attribute values
        overwrite the defaults defined in this class. This mirrors how other pyISSM
        classes support lightweight “inherit defaults then override” workflows.

    Attributes
    ----------
    iscontrol : int
        Flag enabling inversion (0: off, 1: on). When 0, ``check_consistency`` and
        ``marshall_class`` will early-return and write only the minimal fields.
    maxsteps : int
        Maximum number of *gradient evaluations* (outer iterations).
    maxiter : int
        Maximum number of *function evaluations* (forward model runs).
    dxmin : float
        Convergence threshold on the parameter update size. Two iterates closer than
        ``dxmin`` in sup-norm are treated as identical (MATLAB parity).
    dfmin_frac : float
        Expected fractional reduction of the objective during the first step.
        Example: 0.3 corresponds to an expected 30% reduction.
    gttol : float
        Gradient convergence criterion measured relative to the initial gradient:
        ``||g(X)|| / ||g(X0)||``.

    Methods
    -------
    extrude(md)
        No-op. Present for consistency with other classes and to support uniform
        model pipelines that call ``extrude`` on all sub-objects.
    check_consistency(md, solution=None, analyses=None)
        Validate fields and (when possible) verify that ISSM was compiled with M1QN3
        support, which ADM1QN3 depends on.
    marshall_class(fid, prefix, md=None)
        Write this object to the ISSM binary input file using ``execute.WriteData``.

    Notes
    -----
    - ISSM inversion type written by this class is **4** (ADM1QN3).
    - The M1QN3 availability check uses ``tools.wrappers.IssmConfig("_HAVE_M1QN3_")``
      when Python wrappers are available. If wrappers are unavailable, the check is
      skipped with a warning (mirroring patterns elsewhere in pyISSM).
    - This class intentionally does not define cost functions, observation vectors,
      or control parameters: those belong to other inversion formulations in pyISSM.

    Examples
    --------
    Basic setup:

    >>> md.inversion = pyissm.model.classes.inversion.adm1qn3()
    >>> md.inversion.iscontrol = 1
    >>> md.inversion.maxsteps = 30
    >>> md.inversion.gttol = 1e-5

    Inherit settings from another instance:

    >>> base = pyissm.model.classes.inversion.adm1qn3()
    >>> base.maxiter = 60
    >>> md.inversion = pyissm.model.classes.inversion.adm1qn3(other=base)
    """

    def __init__(self, other=None):
        """
        Initialize ADM1QN3 inversion parameters.

        Parameters
        ----------
        other : object, optional
            Object providing attribute overrides. Any attribute on ``other`` with the
            same name as an attribute on this class will replace the default value.
        """
        # Defaults (MATLAB parity)
        self.iscontrol = 0
        self.maxsteps = 20
        self.maxiter = 40
        self.dxmin = 0.1
        self.dfmin_frac = 1.0
        self.gttol = 1e-4

        # Inherit matching fields from provided class/object
        super().__init__(other)

    def __repr__(self):
        """
        Return a detailed, human-readable summary of ADM1QN3 parameters.

        Returns
        -------
        str
            Multi-line representation suitable for console inspection, following the
            same style as other pyISSM classes.
        """
        s = "   adm1qn3 inversion parameters:\n"
        s += f"{class_utils.fielddisplay(self, 'iscontrol', 'is inversion activated?')}\n"
        s += f"{class_utils.fielddisplay(self, 'maxsteps', 'maximum number of iterations (gradient evaluations)')}\n"
        s += f"{class_utils.fielddisplay(self, 'maxiter', 'maximum number of function evaluations (forward runs)')}\n"
        s += f"{class_utils.fielddisplay(self, 'dxmin', 'convergence threshold on parameter updates')}\n"
        s += f"{class_utils.fielddisplay(self, 'dfmin_frac', 'expected reduction of objective during first step')}\n"
        s += f"{class_utils.fielddisplay(self, 'gttol', 'relative gradient convergence criterion')}\n"
        return s

    def __str__(self):
        """
        Return a short identifier string.

        Returns
        -------
        str
            Short class identifier consistent with other inversion classes.
        """
        return "ISSM - inversion.adm1qn3 Class"

    def extrude(self, md):
        """
        Extrude fields to a 3D mesh (no-op).

        ADM1QN3 stores only scalar optimizer settings; there are no node/element fields
        to project when moving from 2D to 3D. This method exists to keep a uniform
        interface across model classes.

        Parameters
        ----------
        md : Model
            ISSM/pyISSM model instance.

        Returns
        -------
        adm1qn3
            Returns ``self`` for chaining.
        """
        return self

    def check_consistency(self, md, solution=None, analyses=None):
        """
        Validate ADM1QN3 fields and installation prerequisites.

        If inversion is disabled (``iscontrol == 0``), this method returns immediately
        without checks. Otherwise it validates:

        - that ISSM was compiled with M1QN3 support (when wrappers are available)
        - scalar bounds and admissible ranges for optimizer parameters

        Parameters
        ----------
        md : Model
            ISSM/pyISSM model instance.
        solution : str, optional
            Present for API compatibility; not used by ADM1QN3.
        analyses : any, optional
            Present for API compatibility; not used by ADM1QN3.

        Returns
        -------
        Model
            The input model (possibly annotated with messages via ``md.check_message``).
        """
        if not self.iscontrol:
            return md

        # Check M1QN3 availability when possible
        if tools.wrappers.check_wrappers_installed():
            if not tools.wrappers.IssmConfig("_HAVE_M1QN3_")[0]:
                md = md.check_message(
                    "M1QN3 has not been installed. ISSM must be reconfigured "
                    "and recompiled with M1QN3."
                )
        else:
            warnings.warn(
                "inversion.adm1qn3.check_consistency: Python wrappers not installed; "
                "skipping M1QN3 availability check."
            )

        class_utils.check_field(md, fieldname="inversion.iscontrol", values=[0, 1])
        class_utils.check_field(md, fieldname="inversion.maxsteps", scalar=True, ge=0)
        class_utils.check_field(md, fieldname="inversion.maxiter", scalar=True, ge=0)
        class_utils.check_field(md, fieldname="inversion.dxmin", scalar=True, gt=0)
        class_utils.check_field(md, fieldname="inversion.dfmin_frac", scalar=True, ge=0.0, le=1.0)
        class_utils.check_field(md, fieldname="inversion.gttol", scalar=True, gt=0)

        return md

    def marshall_class(self, fid, prefix, md=None):
        """
        Serialize ADM1QN3 parameters to ISSM binary input.

        This method writes:

        - ``md.inversion.type = 4`` (ADM1QN3)
        - ``iscontrol`` always
        - if ``iscontrol == 1``, the optimizer parameters:
          ``maxsteps``, ``maxiter``, ``dxmin``, ``dfmin_frac``, ``gttol``

        Parameters
        ----------
        fid : file-like
            Open file handle used for writing ISSM binary input.
        prefix : str
            Prefix used by ISSM for namespacing serialized fields.
        md : Model, optional
            ISSM/pyISSM model instance (unused here, kept for signature consistency).

        Returns
        -------
        None
        """
        # Inversion type ID for ADM1QN3
        execute.WriteData(fid, prefix, name="md.inversion.type", data=4, format="Integer")

        # Always write iscontrol
        execute.WriteData(fid, prefix, obj=self, fieldname="iscontrol", format="Boolean")

        if not self.iscontrol:
            return

        execute.WriteData(fid, prefix, obj=self, fieldname="maxsteps", format="Integer")
        execute.WriteData(fid, prefix, obj=self, fieldname="maxiter", format="Integer")
        execute.WriteData(fid, prefix, obj=self, fieldname="dxmin", format="Double")
        execute.WriteData(fid, prefix, obj=self, fieldname="dfmin_frac", format="Double")
        execute.WriteData(fid, prefix, obj=self, fieldname="gttol", format="Double")
