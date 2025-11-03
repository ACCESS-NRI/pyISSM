"""
Primary class for all ISSM model interactions.
"""
from . import param

class Model():
    """
    ISSM Model Class.

    This class defines a high-level container for all components of an ISSM (Ice Sheet System Model) model.
    It initializes a collection of model components, each of which may store inputs, settings, and results related
    to various aspects of the ice sheet simulation.

    Parameters
    ----------
    None.

    Attributes
    ----------
    mesh : param.mesh.mesh2d()
        Mesh properties.
    mask : param.mask.mask2d()
        Defines grounded and floating elements.
    geometry : param.geometry.geometry2d()
        Surface elevation, bedrock topography, ice thickness, etc.
    constants : param.constants()
        Physical constants.
    smb : param.smb.default()
        Surface mass balance.
    basalforcings : param.basalforcings.default()
        Bed forcings.
    materials : param.materials.ice()
        Material properties.
    damage : param.damage()
        Damage propagation laws.
    friction : param.friction.default()
        Basal friction / drag properties.
    flowequation : param.flowequation()
        Flow equations.
    timestepping : param.timestepping.default()
        Timestepping for transient models.
    initialization : param.initialization()
        Initial guess / state.
    rifts : param.rifts()
        Rifts properties.
    solidearth : param.solidearth.earth()
        Solidearth inputs and settings.
    dsl : param.dsl.default()
        Dynamic sea level.
    debug : param.debug()
        Debugging tools (valgrind, gprof).
    verbose : param.verbose()
        Verbosity level in solve.
    settings : param.issmsettings()
        Settings properties.
    toolkits : None
        PETSc options for each solution.
    cluster : None
        Cluster parameters (number of CPUs, etc.).
    balancethickness : param.balancethickness()
        Parameters for balancethickness solution.
    stressbalance : param.stressbalance()
        Parameters for stressbalance solution.
    groundingline : param.groundingline()
        Parameters for groundingline solution.
    hydrology : param.hydrology.shreve()
        Parameters for hydrology solution.
    masstransport : param.masstransport()
        Parameters for masstransport solution.
    thermal : param.thermal()
        Parameters for thermal solution.
    steadystate : param.steadystate()
        Parameters for steadystate solution.
    transient : param.transient()
        Parameters for transient solution.
    levelset : param.levelset()
        Parameters for moving boundaries (level-set method).
    calving : param.calving.default()
        Parameters for calving.
    frontalforcings : param.frontalforcings.default()
        Parameters for frontalforcings.
    love : param.love.default()
        Parameters for love solution.
    esa : param.esa()
        Parameters for elastic adjustment solution.
    sampling : param.sampling()
        Parameters for stochastic sampler.
    autodiff : param.autodiff()
        Automatic differentiation parameters.
    inversion : param.inversion.default()
        Parameters for inverse methods.
    qmu : param.qmu.default()
        Dakota properties.
    amr : param.amr()
        Adaptive mesh refinement properties.
    results : param.results.default()
        Model results.
    outputdefinition : param.outputdefinition()
        Output definition.
    radaroverlay : param.radaroverlay()
        Radar image for plot overlay.
    miscellaneous : param.miscellaneous()
        Miscellaneous fields.
    stochasticforcing : param.stochasticforcing()
        Stochasticity applied to model forcings.
    """

    def __init__(self):

        ## Initialise all as None
        self.mesh = param.mesh.mesh2d()
        self.mask = param.mask()
        self.geometry = param.geometry()
        self.constants = param.constants()
        self.smb = param.smb.default()
        self.basalforcings = param.basalforcings.default()
        self.materials = param.materials.ice()
        self.damage = param.damage()
        self.friction = param.friction.default()
        self.flowequation = param.flowequation()
        self.timestepping = param.timestepping.default()
        self.initialization = param.initialization()
        self.rifts = param.rifts()
        self.dsl = param.dsl.default()
        self.solidearth = param.solidearth.earth()
        self.debug = param.debug()
        self.verbose = param.verbose()
        self.settings = param.issmsettings()
        self.toolkits = param.toolkits()
        self.cluster = param.cluster.generic()
        self.balancethickness = param.balancethickness()
        self.stressbalance = param.stressbalance()
        self.groundingline = param.groundingline()
        self.hydrology = param.hydrology.shreve()
        self.debris = param.debris()
        self.masstransport = param.masstransport()
        self.thermal = param.thermal()
        self.steadystate = param.steadystate()
        self.transient = param.transient()
        self.levelset = param.levelset()
        self.calving = param.calving.default()
        self.frontalforcings = param.frontalforcings.default()
        self.love = param.love.default()
        self.esa = param.esa()
        self.sampling = param.sampling()
        self.autodiff = param.autodiff()
        self.inversion = param.inversion.default()
        self.qmu = param.qmu.default()
        self.amr = param.amr()
        self.results = param.results.default()
        self.outputdefinition = param.outputdefinition()
        self.radaroverlay = param.radaroverlay()
        self.miscellaneous = param.miscellaneous()
        self.private = param.private()
        self.stochasticforcing = param.stochasticforcing()

    # Define repr
    def __repr__(self):
        # Largely consistent with current MATLAB setup
        s = '%19s %-23s %s' % ('ISSM Model Class', '', '')
        s = '%s\n%s' % (s, '%19s %-23s %s' % ('', '', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('mesh', 'mesh properties', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('mask', 'defines grounded and gloating elements', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('geometry', 'surface elevation, bedrock topography, ice thickness, ...', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('constants', 'physical constants', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('smb', 'surface mass balance', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('basalforcings', 'bed forcings', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('materials', 'material properties', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('damage', 'damage propagation laws', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('friction', 'basal friction / drag properties', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('flowequation', 'flow equations', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('timestepping', 'timestepping for transient models', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('initialization', 'initial guess / state', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('rifts', 'rifts properties', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('solidearth', 'solidearth inputs and settings', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('dsl', 'dynamic sea level', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('debug', 'debugging tools (valgrind, gprof', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('verbose', 'verbosity level in solve', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('settings', 'settings properties', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('toolkits', 'PETSc options for each solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('cluster', 'cluster parameters (number of CPUs...)', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('balancethickness', 'parameters for balancethickness solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('stressbalance', 'parameters for stressbalance solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('groundingline', 'parameters for groundingline solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('hydrology', 'parameters for hydrology solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('masstransport', 'parameters for masstransport solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('thermal', 'parameters fo thermal solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('steadystate', 'parameters for steadystate solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('transient', 'parameters for transient solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('levelset', 'parameters for moving boundaries (level-set method)', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('calving', 'parameters for calving', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('frontalforcings', 'parameters for frontalforcings', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('esa', 'parameters for elastic adjustment solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('sampling', 'parameters for stochastic sampler', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('love', 'parameters for love solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('autodiff', 'automatic differentiation parameters', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('inversion', 'parameters for inverse methods', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('qmu', 'Dakota properties', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('amr', 'adaptive mesh refinement properties', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('outputdefinition', 'output definition', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('results', 'modelresults', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('radaroverlay', 'radar image for plot overlay', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('miscellaneous', 'miscellaneous fields', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('stochasticforcing', 'stochasticity applied to model forcings', ''))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM Model Class'
        return s
    
    def check_message(self, string):
        print('model not consistent: {}'.format(string))
        self.private.isconsistent = False
        return self

    # Define state
    def __getstate__(self):
        return self.__dict__.copy()