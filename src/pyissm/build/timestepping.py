from . import build_utils
from . import class_registry

## ------------------------------------------------------
## timestepping.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    '''
    timestepping.default Class definition
    '''

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
    '''
    timestepping.adaptive Class definition
    '''

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

