import numpy as np
from . import build_utils
from . import class_registry

## ------------------------------------------------------
## inversion.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    '''
    inversion.default Class definition
    '''

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

        s += '{}\n'.format(build_utils.fielddisplay(self, 'iscontrol', 'is inversion activated?'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'incomplete_adjoint', '1: linear viscosity, 0: non - linear viscosity'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'control_parameters', 'ex: {''FrictionCoefficient''}, or {''MaterialsRheologyBbar''}'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'nsteps', 'number of optimization searches'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'cost_functions', 'indicate the type of response for each optimization step'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'cost_functions_coefficients', 'cost_functions_coefficients applied to the misfit of each vertex and for each control_parameter'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'cost_function_threshold', 'misfit convergence criterion. Default is 1%, NaN if not applied'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'maxiter_per_step', 'maximum iterations during each optimization step'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'gradient_scaling', 'scaling factor on gradient direction during optimization, for each optimization step'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'step_threshold', 'decrease threshold for misfit, default is 30%'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'min_parameters', 'absolute minimum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'max_parameters', 'absolute maximum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vx_obs', 'observed velocity x component [m/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vy_obs', 'observed velocity y component [m/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vel_obs', 'observed velocity magnitude [m/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thickness_obs', 'observed thickness [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'surface_obs', 'observed surface elevation [m]'))
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

## ------------------------------------------------------
## inversion.m1qn3
## ------------------------------------------------------
@class_registry.register_class
class m1qn3(class_registry.manage_state):
    '''
    inversion.m1qn3 Class definition
    '''

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

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   m1qn3 inversion parameters:\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'iscontrol', 'is inversion activated?'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'incomplete_adjoint', '1: linear viscosity, 0: non - linear viscosity'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'control_parameters', 'ex: [''FrictionCoefficient''], or [''MaterialsRheologyBbar'']'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'control_scaling_factors', 'order of magnitude of each control (useful for multi - parameter optimization)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'maxsteps', 'maximum number of iterations (gradient computation)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'maxiter', 'maximum number of Function evaluation (forward run)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'dxmin', 'convergence criterion: two points less than dxmin from eachother (sup - norm) are considered identical'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'dfmin_frac', 'expected reduction of J during the first step (e.g., 0.3=30% reduction in cost function)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'gttol', '||g(X)||/||g(X0)|| (g(X0): gradient at initial guess X0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'cost_functions', 'indicate the type of response for each optimization step'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'cost_functions_coefficients', 'cost_functions_coefficients applied to the misfit of each vertex and for each control_parameter'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'min_parameters', 'absolute minimum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'max_parameters', 'absolute maximum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vx_obs', 'observed velocity x component [m / yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vy_obs', 'observed velocity y component [m / yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vel_obs', 'observed velocity magnitude [m / yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thickness_obs', 'observed thickness [m]'))
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

## ------------------------------------------------------
## inversion.tao
## ------------------------------------------------------
@class_registry.register_class
class tao(class_registry.manage_state):
    '''
    inversion.tao Class definition
    '''

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
        self.algorithm = 'blmvm' # TODO: Add check for PETSC version from IssmConfig
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

        s += '{}\n'.format(build_utils.fielddisplay(self, 'iscontrol', 'is inversion activated?'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'incomplete_adjoint', '1: linear viscosity, 0: non - linear viscosity'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'control_parameters', 'ex: {''FrictionCoefficient''}, or {''MaterialsRheologyBbar''}'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'maxsteps', 'maximum number of iterations (gradient computation)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'maxiter', 'maximum number of Function evaluation (forward run)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'fatol', 'convergence criterion: f(X) - f(X * ) (X: current iteration, X * : "true" solution, f: cost function)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'frtol', 'convergence criterion: |f(X) - f(X * )| / |f(X * )|'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'gatol', 'convergence criterion: ||g(X)|| (g: gradient of the cost function)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'grtol', 'convergence criterion: ||g(X)|| / |f(X)|'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'gttol', 'convergence criterion: ||g(X)|| / ||g(X0)|| (g(X0): gradient at initial guess X0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'algorithm', 'minimization algorithm: ''tao_blmvm'', ''tao_cg'', ''tao_lmvm'''))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'cost_functions', 'indicate the type of response for each optimization step'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'cost_functions_coefficients', 'cost_functions_coefficients applied to the misfit of each vertex and for each control_parameter'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'min_parameters', 'absolute minimum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'max_parameters', 'absolute maximum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vx_obs', 'observed velocity x component [m / yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vy_obs', 'observed velocity y component [m / yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vel_obs', 'observed velocity magnitude [m / yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thickness_obs', 'observed thickness [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'surface_obs', 'observed surface elevation [m]'))
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