"""
Primary class for all ISSM model interactions.
"""
from . import build

class Model():
    """
    ISSM Model Class.

    This class defines a high-level container for all components of an ISSM (Ice Sheet System Model) model.
    It initializes a collection of model components, each of which may store inputs, settings, and results related
    to various aspects of the ice sheet simulation.

    Parameters
    ----------
    args : optional
        Optional arguments passed to the model.

    Attributes
    ----------
    mesh : SimpleNamespace or None
        Mesh properties.
    mask : SimpleNamespace or None
        Defines grounded and floating elements.
    geometry : SimpleNamespace or None
        Surface elevation, bedrock topography, ice thickness, etc.
    constants : SimpleNamespace or None
        Physical constants.
    smb : SimpleNamespace or None
        Surface mass balance.
    basalforcings : SimpleNamespace or None
        Bed forcings.
    materials : SimpleNamespace or None
        Material properties.
    damage : SimpleNamespace or None
        Damage propagation laws.
    friction : SimpleNamespace or None
        Basal friction / drag properties.
    flowequation : SimpleNamespace or None
        Flow equations.
    timestepping : SimpleNamespace or None
        Timestepping for transient models.
    initialization : SimpleNamespace or None
        Initial guess / state.
    rifts : SimpleNamespace or None
        Rifts properties.
    solidearth : SimpleNamespace or None
        Solidearth inputs and settings.
    dsl : SimpleNamespace or None
        Dynamic sea level.
    debug : SimpleNamespace or None
        Debugging tools (valgrind, gprof).
    verbose : SimpleNamespace or None
        Verbosity level in solve.
    settings : SimpleNamespace or None
        Settings properties.
    toolkits : SimpleNamespace or None
        PETSc options for each solution.
    cluster : SimpleNamespace or None
        Cluster parameters (number of CPUs, etc.).
    balancethickness : SimpleNamespace or None
        Parameters for balancethickness solution.
    stressbalance : SimpleNamespace or None
        Parameters for stressbalance solution.
    groundingline : SimpleNamespace or None
        Parameters for groundingline solution.
    hydrology : SimpleNamespace or None
        Parameters for hydrology solution.
    masstransport : SimpleNamespace or None
        Parameters for masstransport solution.
    thermal : SimpleNamespace or None
        Parameters for thermal solution.
    steadystate : SimpleNamespace or None
        Parameters for steadystate solution.
    transient : SimpleNamespace or None
        Parameters for transient solution.
    levelset : SimpleNamespace or None
        Parameters for moving boundaries (level-set method).
    calving : SimpleNamespace or None
        Parameters for calving.
    frontalforcings : SimpleNamespace or None
        Parameters for frontalforcings.
    esa : SimpleNamespace or None
        Parameters for elastic adjustment solution.
    sampling : SimpleNamespace or None
        Parameters for stochastic sampler.
    love : SimpleNamespace or None
        Parameters for love solution.
    autodiff : SimpleNamespace or None
        Automatic differentiation parameters.
    inversion : SimpleNamespace or None
        Parameters for inverse methods.
    qmu : SimpleNamespace or None
        Dakota properties.
    amr : SimpleNamespace or None
        Adaptive mesh refinement properties.
    results : SimpleNamespace or None
        Model results.
    outputdefinition : SimpleNamespace or None
        Output definition.
    radaroverlay : SimpleNamespace or None
        Radar image for plot overlay.
    miscellaneous : SimpleNamespace or None
        Miscellaneous fields.
    stochasticforcing : SimpleNamespace or None
        Stochasticity applied to model forcings.
    """

    def __init__(self, args=None):

        self.args = args

        ## Initialise all as None
        self.mesh = build.mesh.mesh2d()
        self.mask = build.mask()
        self.geometry = build.geometry()
        self.constants = build.constants()
        self.smb = build.smb.default()
        self.basalforcings = build.basalforcings.default()
        self.materials = build.materials.ice()
        self.damage = build.damage()
        self.friction = build.friction.default()
        self.flowequation = build.flowequation()
        self.timestepping = build.timestepping.default()
        self.initialization = build.initialization()
        self.rifts = build.rifts()
        self.dsl = build.dsl.default()
        self.solidearth = build.solidearth.earth()
        self.debug = build.debug()
        self.verbose = None
        self.settings = build.issmsettings()
        self.toolkits = None
        self.cluster = None
        self.balancethickness = build.balancethickness()
        self.stressbalance = build.stressbalance()
        self.groundingline = build.groundingline()
        self.hydrology = build.hydrology.shreve()
        self.debris = build.debris()
        self.masstransport = build.masstransport()
        self.thermal = build.thermal()
        self.steadystate = build.steadystate()
        self.transient = build.transient()
        self.levelset = build.levelset()
        self.calving = build.calving.default()
        self.frontalforcings = build.frontalforcings.default()
        self.love = build.love.default()
        self.esa = build.esa()
        self.sampling = build.sampling()
        self.autodiff = build.autodiff()
        self.inversion = build.inversion.default()
        self.qmu = build.qmu.default()
        self.amr = build.amr()
        self.results = build.results.default()
        self.outputdefinition = build.outputdefinition()
        self.radaroverlay = build.radaroverlay()
        self.miscellaneous = build.miscellaneous()
        self.private = build.private()
        self.stochasticforcing = build.stochasticforcing()

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