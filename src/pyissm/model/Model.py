"""
Primary class for all ISSM model interactions.
"""
from pyissm.model import classes

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
    mesh : classes.mesh.mesh2d()
        Mesh properties.
    mask : classes.mask.mask2d()
        Defines grounded and floating elements.
    geometry : classes.geometry.geometry2d()
        Surface elevation, bedrock topography, ice thickness, etc.
    constants : classes.constants()
        Physical constants.
    smb : classes.smb.default()
        Surface mass balance.
    basalforcings : classes.basalforcings.default()
        Bed forcings.
    materials : classes.materials.ice()
        Material properties.
    damage : classes.damage()
        Damage propagation laws.
    friction : classes.friction.default()
        Basal friction / drag properties.
    flowequation : classes.flowequation()
        Flow equations.
    timestepping : classes.timestepping.default()
        Timestepping for transient models.
    initialization : classes.initialization()
        Initial guess / state.
    rifts : classes.rifts()
        Rifts properties.
    solidearth : classes.solidearth.earth()
        Solidearth inputs and settings.
    dsl : classes.dsl.default()
        Dynamic sea level.
    debug : classes.debug()
        Debugging tools (valgrind, gprof).
    verbose : classes.verbose()
        Verbosity level in solve.
    settings : classes.issmsettings()
        Settings properties.
    toolkits : None
        PETSc options for each solution.
    cluster : None
        Cluster parameters (number of CPUs, etc.).
    balancethickness : classes.balancethickness()
        Parameters for balancethickness solution.
    stressbalance : classes.stressbalance()
        Parameters for stressbalance solution.
    groundingline : classes.groundingline()
        Parameters for groundingline solution.
    hydrology : classes.hydrology.shreve()
        Parameters for hydrology solution.
    masstransport : classes.masstransport()
        Parameters for masstransport solution.
    thermal : classes.thermal()
        Parameters for thermal solution.
    steadystate : classes.steadystate()
        Parameters for steadystate solution.
    transient : classes.transient()
        Parameters for transient solution.
    levelset : classes.levelset()
        Parameters for moving boundaries (level-set method).
    calving : classes.calving.default()
        Parameters for calving.
    frontalforcings : classes.frontalforcings.default()
        Parameters for frontalforcings.
    love : classes.love.default()
        Parameters for love solution.
    esa : classes.esa()
        Parameters for elastic adjustment solution.
    sampling : classes.sampling()
        Parameters for stochastic sampler.
    autodiff : classes.autodiff()
        Automatic differentiation parameters.
    inversion : classes.inversion.default()
        Parameters for inverse methods.
    qmu : classes.qmu.default()
        Dakota properties.
    amr : classes.amr()
        Adaptive mesh refinement properties.
    results : classes.results.default()
        Model results.
    outputdefinition : classes.outputdefinition()
        Output definition.
    radaroverlay : classes.radaroverlay()
        Radar image for plot overlay.
    miscellaneous : classes.miscellaneous()
        Miscellaneous fields.
    stochasticforcing : classes.stochasticforcing()
        Stochasticity applied to model forcings.
    """

    def __init__(self):

        ## Initialise all as None
        self.mesh = classes.mesh.mesh2d()
        self.mask = classes.mask()
        self.geometry = classes.geometry()
        self.constants = classes.constants()
        self.smb = classes.smb.default()
        self.basalforcings = classes.basalforcings.default()
        self.materials = classes.materials.ice()
        self.damage = classes.damage()
        self.friction = classes.friction.default()
        self.flowequation = classes.flowequation()
        self.timestepping = classes.timestepping.default()
        self.initialization = classes.initialization()
        self.rifts = classes.rifts()
        self.dsl = classes.dsl.default()
        self.solidearth = classes.solidearth.earth()
        self.debug = classes.debug()
        self.verbose = classes.verbose()
        self.settings = classes.issmsettings()
        self.toolkits = classes.toolkits()
        self.cluster = classes.cluster.generic()
        self.balancethickness = classes.balancethickness()
        self.stressbalance = classes.stressbalance()
        self.groundingline = classes.groundingline()
        self.hydrology = classes.hydrology.shreve()
        self.debris = classes.debris()
        self.masstransport = classes.masstransport()
        self.thermal = classes.thermal()
        self.steadystate = classes.steadystate()
        self.transient = classes.transient()
        self.levelset = classes.levelset()
        self.calving = classes.calving.default()
        self.frontalforcings = classes.frontalforcings.default()
        self.love = classes.love.default()
        self.esa = classes.esa()
        self.sampling = classes.sampling()
        self.autodiff = classes.autodiff()
        self.inversion = classes.inversion.default()
        self.qmu = classes.qmu.default()
        self.amr = classes.amr()
        self.results = classes.results.default()
        self.outputdefinition = classes.outputdefinition()
        self.radaroverlay = classes.radaroverlay()
        self.miscellaneous = classes.miscellaneous()
        self.private = classes.private()
        self.stochasticforcing = classes.stochasticforcing()

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
        """
        Notify about a model consistency error, update internal state, and return the instance.

        This method prints a formatted consistency error message to standard output,
        marks the instance as inconsistent by setting ``self.private.isconsistent``
        to ``False``, and returns the instance to allow for method chaining.

        Parameters
        ----------
        string : str
            Human-readable description of the consistency error. This will be inserted
            into the printed message: ``Model consistency error: {string}``.

        Returns
        -------
        self
            The same instance on which the method was called, enabling fluent/chained
            calls.

        Notes
        -----
        This method has the side effect of mutating the instance state (``self.private.isconsistent``),
        and it performs output via ``print``. It does not raise exceptions.

        Examples
        --------
        >>> obj.check_message("missing parameter")
        Model consistency error: missing parameter
        >>> obj.private.isconsistent
        False
        """
        print(f'Model consistency error: {string}')
        self.private.isconsistent = False
        return self
    
    def model_class_names(self):
        """
        Return a sorted list of registered model class attribute names.

        The method inspects the instance attributes and returns those whose
        classes are registered in ``classes.class_registry.CLASS_REGISTRY``.

        Returns
        -------
        list of str
            Sorted list of attribute names corresponding to registered model classes.
        """
        registered_classes = set(classes.class_registry.CLASS_REGISTRY.values())
        names = [
            name for name, obj in vars(self).items()
            if obj.__class__ in registered_classes
        ]
        return sorted(names)

    # Define state
    def __getstate__(self):
        return self.__dict__.copy()