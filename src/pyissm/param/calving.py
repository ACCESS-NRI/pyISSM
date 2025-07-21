import numpy as np
from . import param_utils
from . import class_registry

## ------------------------------------------------------
## calving.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    """
    Default calving parameters class for ISSM.

    This class encapsulates the default parameters for calving in the ISSM (Ice Sheet System Model) framework.
    It defines the calving rate parameter.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    calvingrate : ndarray, default=np.nan
        Calving rate at given location [m/a].

    Methods
    -------
    __init__(self, other=None)
        Initializes the calving parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the calving parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.calving = pyissm.build.calving.default()
    md.calving.calvingrate = np.zeros((md.mesh.numberofvertices,))
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.calvingrate = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'calvingrate', 'calving rate at given location [m/a]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.default Class'
        return s

## ------------------------------------------------------
## calving.crevassedepth
## ------------------------------------------------------
@class_registry.register_class
class crevassedepth(class_registry.manage_state):
    """
    Crevasse depth calving parameters class for ISSM.

    This class encapsulates the parameters for the crevasse depth calving model in the ISSM (Ice Sheet System Model) framework.
    It defines parameters related to crevasse opening stress, threshold, and water height.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    crevasse_opening_stress : int, default=1
        0: stress only in the ice-flow direction, 1: maximum principal stress.
    crevasse_threshold : float, default=1.0
        Ratio of full thickness to calve (e.g. 0.75 is for 75% of the total ice thickness).
    water_height : float, default=0.0
        Water height in the crevasse [m].

    Methods
    -------
    __init__(self, other=None)
        Initializes the crevasse depth calving parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the crevasse depth calving parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.calving = pyissm.build.calving.crevassedepth()
    md.calving.crevasse_opening_stress = 1.0
    md.calving.crevasse_threshold = 0.75
    md.calving.water_height = 10.0
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.crevasse_opening_stress = 1
        self.crevasse_threshold = 1.
        self.water_height = 0.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving Pi parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'crevasse_opening_stress', '0: stress only in the ice-flow direction, 1: max principal'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'crevasse_threshold', 'ratio of full thickness to calve (e.g. 0.75 is for 75% of the total ice thickness)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'water_height', 'water height in the crevasse [m]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.crevassedepth Class'
        return s

## ------------------------------------------------------
## calving.dev
## ------------------------------------------------------
@class_registry.register_class
class dev(class_registry.manage_state):
    """
    Development calving parameters class for ISSM.

    This class encapsulates the parameters for the development calving model in the ISSM (Ice Sheet System Model) framework.
    It defines stress thresholds for grounded and floating ice.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    stress_threshold_groundedice : float, default=1e6
        Maximum stress threshold applied to grounded ice [Pa].
    stress_threshold_floatingice : float, default=150e3
        Maximum stress threshold applied to floating ice [Pa].

    Methods
    -------
    __init__(self, other=None)
        Initializes the development calving parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the development calving parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.calving = pyissm.build.calving.dev()
    md.calving.stress_threshold_groundedice = 2e6
    md.calving.stress_threshold_floatingice = 200e3
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.stress_threshold_groundedice = 1e6
        self.stress_threshold_floatingice = 150e3

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving Pi parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'stress_threshold_groundedice', 'sigma_max applied to grounded ice only [Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'stress_threshold_floatingice', 'sigma_max applied to floating ice only [Pa]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.dev Class'
        return s

## ------------------------------------------------------
## calving.levermann
## ------------------------------------------------------
@class_registry.register_class
class levermann(class_registry.manage_state):
    """
    Levermann calving parameters class for ISSM.

    This class encapsulates the parameters for the Levermann calving model in the ISSM (Ice Sheet System Model) framework.
    It defines the proportionality coefficient used in the Levermann calving law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    coeff : float, default=2e13
        Proportionality coefficient in the Levermann calving model.

    Methods
    -------
    __init__(self, other=None)
        Initializes the Levermann calving parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the Levermann calving parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.calving = pyissm.build.calving.levermann()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.coeff = 2e13

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving Levermann parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'coeff', 'proportionality coefficient in Levermann model'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.levermann Class'
        return s

## ------------------------------------------------------
## calving.minthickness
## ------------------------------------------------------
@class_registry.register_class
class minthickness(class_registry.manage_state):
    """
    Minimum thickness calving parameters class for ISSM.

    This class encapsulates the parameters for the minimum thickness calving model in the ISSM (Ice Sheet System Model) framework.
    It defines the minimum ice thickness below which no ice is allowed.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    min_thickness : float, default=100.
        Minimum thickness below which no ice is allowed [m].

    Methods
    -------
    __init__(self, other=None)
        Initializes the minimum thickness calving parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the minimum thickness calving parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.calving = pyissm.build.calving.minthickness()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.min_thickness = 100.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving Minimum thickness:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'min_thickness', 'minimum thickness below which no ice is allowed'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.minthickness Class'
        return s

## ------------------------------------------------------
## calving.parameterization
## ------------------------------------------------------
@class_registry.register_class
class parameterization(class_registry.manage_state):
    """
    Parameterization calving parameters class for ISSM.

    This class encapsulates the parameters for the parameterization calving model in the ISSM (Ice Sheet System Model) framework.
    It defines parameters controlling the calving rate as a function of ice thickness, velocity, and other factors.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    min_thickness : float, default=0.
        Minimum thickness below which no ice is allowed [m].
    use_param : float, default=0.
        Parameterization mode selector:
            -1: use frontal ablation rate,
             0: linear function,
             1: tanh function,
             2: tanh(thickness),
             3: tanh(normalized velocity),
             4: tanh(truncated velocity),
             5: linear(truncated velocity).
    theta : float, default=0.
        Amplifier parameter for the calving law.
    alpha : float, default=0.
        Slope parameter for the calving law.
    xoffset : float, default=0.
        Offset in x-axis for the calving law.
    yoffset : float, default=0.
        Offset in y-axis for the calving law.
    vel_upperbound : float, default=6000.
        Upper bound of ice velocity to reduce the calving rate [m/a].
    vel_threshold : float, default=0.
        Threshold of ice velocity to reduce the calving rate [m/a].
    vel_lowerbound : float, default=0.
        Lower bound of ice velocity to reduce the calving rate [m/a].

    Methods
    -------
    __init__(self, other=None)
        Initializes the parameterization calving parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the parameterization calving parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.calving = pyissm.build.calving.parameterization()
    md.calving.min_thickness = 50.
    md.calving.use_param = 1
    md.calving.theta = 2.0
    md.calving.alpha = 0.01
    md.calving.xoffset = 10.0
    md.calving.yoffset = 0.5
    md.calving.vel_upperbound = 5000.
    md.calving.vel_threshold = 1000.
    md.calving.vel_lowerbound = 100.
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.min_thickness = 0.
        self.use_param = 0.
        self.theta = 0.
        self.alpha = 0.
        self.xoffset = 0.
        self.yoffset = 0.
        self.vel_upperbound = 6000.
        self.vel_threshold = 0.
        self.vel_lowerbound = 0.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving test parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'min_thickness', 'minimum thickness below which no ice is allowed [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'use_param', '-1 - just use frontal ablation rate, 0 - f(x) = y_{o} + \alpha (x+x_{o}), 1 - f(x)=y_{o}-\frac{\theta}{2}\tanh(\alpha(x+x_{o})), 2 - tanh(thickness), 3 - tanh(normalized vel), 4 - tanh(truncated vel), 5 - linear(truncated vel)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'theta', 'the amplifier'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'alpha', 'the slope'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'xoffset', 'offset in x-axis'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'yoffset', 'offset in y-axis'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vel_lowerbound', 'lowerbound of ice velocity to reduce the calving rate [m/a]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vel_threshold', 'threshold of ice velocity to reduce the calving rate [m/a]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vel_upperbound', 'upperbound of ice velocity to reduce the calving rate [m/a]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.parameterization Class'
        return s

## ------------------------------------------------------
## calving.vonmises
## ------------------------------------------------------
@class_registry.register_class
class vonmises(class_registry.manage_state):
    """
    Von Mises calving parameters class for ISSM.

    This class encapsulates the parameters for the Von Mises calving model in the ISSM (Ice Sheet System Model) framework.
    It defines stress thresholds for grounded and floating ice, as well as the minimum ice thickness below which no ice is allowed.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    stress_threshold_groundedice : float, default=0
        Maximum Von Mises stress threshold applied to grounded ice [Pa].
    stress_threshold_floatingice : float, default=0
        Maximum Von Mises stress threshold applied to floating ice [Pa].
    min_thickness : float, default=0.
        Minimum thickness below which no ice is allowed [m].

    Methods
    -------
    __init__(self, other=None)
        Initializes the Von Mises calving parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the Von Mises calving parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.calving = pyissm.build.calving.vonmises()
    md.calving.stress_threshold_groundedice = 1e6
    md.calving.stress_threshold_floatingice = 150e3
    md.calving.min_thickness = 50.
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.stress_threshold_groundedice = 0
        self.stress_threshold_floatingice = 0
        self.min_thickness = 0.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving VonMises parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'stress_threshold_groundedice', 'sigma_max applied to grounded ice only [Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'stress_threshold_floatingice', 'sigma_max applied to floating ice only [Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'min_thickness', 'minimum thickness below which no ice is allowed [m]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.vonmises Class'
        return s