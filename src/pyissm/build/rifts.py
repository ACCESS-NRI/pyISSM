from . import build_utils
from . import class_registry

@class_registry.register_class
class rifts(class_registry.manage_state):
    '''
    rifts Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.riftstruct = 'Rift structure'
        self.riftproperties = 'Rift properties'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   rift parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'riftstruct', 'structure containing all rift information (vertices coordinates, segments, type of melange, ...)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'riftproperties', 'rift properties'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - rifts Class'
        return s

