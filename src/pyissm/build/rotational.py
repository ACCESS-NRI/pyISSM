from . import build_utils
from . import class_registry

@class_registry.register_class
class rotational(class_registry.manage_state):
    '''
    rotational Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.equatorialmoi = 0
        self.polarmoi = 0
        self.langularvelocity = 0


        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   rotational parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'equatorialmoi', 'mean equatorial moment of inertia [kg m^2]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'polarmoi', 'polar moment of inertia [kg m^2]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'angularvelocity', 'mean rotational velocity of earth [per second]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - rotational Class'
        return s

