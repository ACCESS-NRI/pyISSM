from . import build_utils
from . import class_registry

@class_registry.register_class
class debug(class_registry.manage_state):
    '''
    debug Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.valgrind = 0
        self.gprof = 0
        self.profiling = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   debug parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'valgrind', 'use valgrind to debug (0 or 1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'gprof', 'use gnu - profiler to find out where the time is spent'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'profiling', 'enables profiling (memory, flops, time)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - debug Class'
        return s

