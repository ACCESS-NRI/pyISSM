from . import build_utils
from . import class_registry

## ------------------------------------------------------
## timestepping.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    """
    Default fixed time stepping configuration for ISSM transient simulations.

    This class configures fixed time step parameters for transient ice sheet simulations.
    It provides control over simulation duration, time step size, and forcing interpolation
    methods for consistent temporal evolution of the ice sheet system.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    start_time : float, default=0
        Simulation starting time [yr].
    final_time : float, default=5
        Final time to stop the simulation [yr].
    time_step : float, default=0.5
        Length of time steps [yr].
    interp_forcing : int, default=1
        Interpolate in time between requested forcing values (0 or 1).
    average_forcing : int, default=0
        Average in time if there are several forcing values between steps (0 or 1).
    cycle_forcing : int, default=0
        Cycle through forcing data (0 or 1).
    coupling_time : float, default=0
        Length of coupling time steps with ocean model [yr].

    Methods
    -------
    __init__(self, other=None)
        Initializes the fixed timestepping parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the timestepping parameters.
    __str__(self)
        Returns a short string identifying the class.

    Notes
    -----
    Fixed time stepping uses constant time step intervals throughout the simulation.
    This approach provides predictable temporal resolution but may not be optimal
    for all simulation conditions.

    Forcing interpolation options:
    - interp_forcing=1: Linear interpolation between forcing time points
    - average_forcing=1: Time averaging when multiple forcing values exist per step
    - cycle_forcing=1: Repeat forcing data cyclically

    Examples
    --------
    md.timestepping = pyissm.build.timestepping.default()
    md.timestepping.start_time = 0
    md.timestepping.final_time = 100
    md.timestepping.time_step = 1.0
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.start_time = 0
        self.final_time = 5
        self.time_step = 0.5
        self.interp_forcing = 1
        self.average_forcing = 0
        self.cycle_forcing = 0
        self.coupling_time = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   timestepping parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'start_time', 'simulation starting time [yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'final_time', 'final time to stop the simulation [yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'time_step', 'length of time steps [yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'interp_forcing', 'interpolate in time between requested forcing values? (0 or 1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'average_forcing', 'average in time if there are several forcing values between steps? (0 or 1, default is 0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'cycle_forcing', 'cycle through forcing? (0 or 1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'coupling_time', 'length of coupling time steps with ocean model [yr]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - timestepping.default Class'
        return s

## ------------------------------------------------------
## timestepping.adaptive
## ------------------------------------------------------
@class_registry.register_class
class adaptive(class_registry.manage_state):
    """
    Adaptive time stepping configuration for ISSM transient simulations.

    This class configures adaptive time step parameters for transient ice sheet simulations.
    The time step size is automatically adjusted based on CFL (Courant-Friedrichs-Lewy)
    conditions to maintain numerical stability while optimizing computational efficiency.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    start_time : float, default=0.0
        Simulation starting time [yr].
    final_time : float, default=100.0
        Final time to stop the simulation [yr].
    time_step_min : float, default=0.01
        Minimum length of time steps [yr].
    time_step_max : float, default=10.0
        Maximum length of time steps [yr].
    cfl_coefficient : float, default=0.5
        Coefficient applied to CFL condition for time step calculation.
    interp_forcing : int, default=1
        Interpolate in time between requested forcing values (0 or 1).
    average_forcing : int, default=0
        Average in time if there are several forcing values between steps (0 or 1).
    cycle_forcing : int, default=0
        Cycle through forcing data (0 or 1).
    coupling_time : float, default=0.0
        Coupling time steps with ocean model [yr].

    Methods
    -------
    __init__(self, other=None)
        Initializes the adaptive timestepping parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the adaptive timestepping parameters.
    __str__(self)
        Returns a short string identifying the class.

    Notes
    -----
    Adaptive time stepping automatically adjusts time step size based on the CFL condition:
    Δt ≤ CFL_coefficient × (grid_spacing / max_velocity)
    
    This approach provides:
    - Numerical stability by respecting CFL limits
    - Computational efficiency by using larger steps when possible
    - Automatic handling of varying flow conditions

    The time step is constrained between time_step_min and time_step_max to prevent
    excessively small or large time steps that could cause numerical issues.

    CFL coefficient recommendations:
    - 0.5: Conservative, very stable (default)
    - 0.8-0.9: More aggressive, faster computation
    - >1.0: May cause instability

    Examples
    --------
    md.timestepping = pyissm.build.timestepping.adaptive()
    md.timestepping.time_step_min = 0.001
    md.timestepping.time_step_max = 5.0
    md.timestepping.cfl_coefficient = 0.8
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.start_time = 0.
        self.final_time = 100.
        self.time_step_min = 0.01
        self.time_step_max = 10.
        self.cfl_coefficient = 0.5
        self.interp_forcing = 1
        self.average_forcing = 0
        self.cycle_forcing = 0
        self.coupling_time = 0.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   timestepping.adaptive parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'start_time', 'simulation starting time [yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, "start_time", "simulation starting time [yr]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "final_time", "final time to stop the simulation [yr]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "time_step_min", "minimum length of time steps [yr]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "time_step_max", "maximum length of time steps [yr]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "cfl_coefficient", "coefficient applied to cfl condition"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "interp_forcing", "interpolate in time between requested forcing values ? (0 or 1)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'average_forcing', 'average in time if there are several forcing values between steps? (0 or 1, default is 0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, "cycle_forcing", "cycle through forcing ? (0 or 1)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "coupling_time", "coupling time steps with ocean model [yr]"))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - timestepping.adaptive Class'
        return s

