from . import build_utils
from . import class_registry

@class_registry.register_class
class steadystate(class_registry.manage_state):
    '''
    steadystate Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.reltol = 0.01
        self.maxiter = 100
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   steadystate solution parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'reltol', 'relative tolerance criterion'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'maxiter', 'maximum number of iterations'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional requested outputs'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - steadystate Class'
        return s

