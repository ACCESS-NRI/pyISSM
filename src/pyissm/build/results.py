from . import build_utils
from . import class_registry

## ------------------------------------------------------
## results.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    """
    Default results container class for ISSM.

    This class serves as the default results container in the ISSM (Ice Sheet System Model) framework.
    It stores simulation results and provides methods for displaying and accessing the stored data.
    The class dynamically stores results as attributes and provides formatted string representations
    of the results structure.

    Parameters
    ----------
    None

    Methods
    -------
    __init__(self)
        Initializes the default results container.
    __repr__(self)
        Returns a detailed string representation of the results structure showing field names and sizes.
    __str__(self)
        Returns a short string identifying the class.

    Notes
    -----
    This class dynamically stores results as attributes. The actual attributes depend on the 
    simulation type and requested outputs. Common attributes may include velocity fields, 
    thickness, pressure, temperature, and other solution variables.

    Examples
    --------
    results = pyissm.build.results.default()
    """

    # Initialise with default parameters
    def __init__(self):
        pass

    # Define repr
    def __repr__(self):  #{{{
        s = ''
        for key, value in self.__dict__.items():
            # TODO: Is this check necessary? resultsdakota is a separate class now.
            if isinstance(value, resultsdakota):
                lengthvalue = 1
            else:
                try:
                    lengthvalue = len(value)
                except TypeError:
                    lengthvalue = 1
            s += '    {}: [1x{} struct]\n'.format(key, lengthvalue)
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - results.default Class'
        return s

## ------------------------------------------------------
## results.resultsdakota
## ------------------------------------------------------
@class_registry.register_class
class resultsdakota(class_registry.manage_state):
    """
    Results container class for Dakota-based ISSM runs.

    This class is designed to store and manage results from Dakota uncertainty quantification or optimization
    runs within the ISSM (Ice Sheet System Model) framework. It dynamically stores results as attributes,
    which may include lists of results for each Dakota evaluation, summary statistics, or other relevant data.

    Parameters
    ----------
    None

    Methods
    -------
    __init__(self)
        Initializes the resultsdakota container.
    __repr__(self)
        Returns a string representation of the results structure showing field names and summary information.
    __str__(self)
        Returns a short string identifying the class.

    Notes
    -----
    This class is typically used when ISSM is run in conjunction with Dakota for parameter studies,
    uncertainty quantification, or optimization. The actual attributes depend on the Dakota study
    configuration and requested outputs.

    Examples
    --------
    results = pyissm.build.results.resultsdakota()
    """

    # Initialise with default parameters
    def __init__(self):
        pass

    # Define repr
    def __repr__(self):
        s = ''
        for key, value in self.__dict__.items():
            s += '    {}: '.format(key)
            if isinstance(value, list):
                s += '[{} element list]'.format(len(value))
            else:
                s += '{}'.format(value)
            s += '\n'
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - resultsdakota Class'
        return s

## ------------------------------------------------------
## results.solution
## ------------------------------------------------------
@class_registry.register_class
class solution(class_registry.manage_state):
    """
    Results container class for ISSM solution steps.

    This class is designed to store and manage the results of solution steps within the ISSM (Ice Sheet System Model)
    framework. Each instance contains a list of solution steps, where each step holds the results for a particular
    time or iteration in the simulation.

    Parameters
    ----------
    args : list, optional
        If provided, should be a list of solutionstep instances. If not provided, initializes with a single default solutionstep.

    Attributes
    ----------
    steps : list of solutionstep
        List containing the solution steps for the simulation.

    Methods
    -------
    __init__(self, *args)
        Initializes the solution container, optionally with a list of solutionstep instances.
    __repr__(self)
        Returns a string representation of the solution structure, showing field names and values.
    __str__(self)
        Returns a short string identifying the class.

    Notes
    -----
    This class is typically used to organize results from time-dependent or iterative ISSM simulations.
    Each solutionstep instance in the steps list contains the results for a single step.

    Examples
    --------
    results = pyissm.build.results.solution()
    """

    # Initialise with default parameters
    def __init__(self, *args):
        self.steps = None
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, list):
                self.steps = arg
            else:
                raise Exception('solution class error: if initializing with an argument, that argument should be an empty list or a list of instances of solutionstep')
        else:
            self.steps = [solutionstep()]

    # Define repr
    def __repr__(self):  #{{{
        s = ''
        numsteps = len(self.steps)
        if numsteps == 1:
            for key, value in self.steps[0].__dict__.items():
                s += '    {}: {}\n'.format(key, value)
        else:
            s = '  1x{} struct array with fields:\n'.format(numsteps)
            s += '\n'
            for fieldname in self.steps[0].getfieldnames():
                s += '    {}\n'.format(fieldname)
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - solution Class'
        return s

## ------------------------------------------------------
## results.solutionstep
## ------------------------------------------------------
@class_registry.register_class
class solutionstep(class_registry.manage_state):
    """
    Results container class for a single ISSM solution step.

    This class is designed to store and manage the results for a single solution step within the ISSM (Ice Sheet System Model)
    framework. Each instance holds the results for a particular time or iteration in the simulation, such as velocity, thickness,
    temperature, or other relevant fields.

    Parameters
    ----------
    None
    
    Attributes
    ----------
    (Dynamic)
        Attributes are dynamically assigned based on the simulation outputs for this step. Typical attributes may include
        velocity, thickness, pressure, temperature, etc.

    Methods
    -------
    __init__(self)
        Initializes the solutionstep container.
    __repr__(self)
        Returns a string representation of the solutionstep structure, showing field names and values.
    __str__(self)
        Returns a short string identifying the class.

    Notes
    -----
    This class is typically used as an element of the steps list in the solution class, representing a single time step or iteration.

    Examples
    --------
    step = pyissm.build.results.solutionstep()
    """

    # Initialise with default parameters
    def __init__(self):
        pass

    # Define repr
    def __repr__(self):
        s = ''
        width = build_utils.getlongestfieldname(self)
        for key, value in self.__dict__.items():
            s += '    {:{width}s}: {}\n'.format(key, value, width=width)
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - solutionstep Class'
        return s