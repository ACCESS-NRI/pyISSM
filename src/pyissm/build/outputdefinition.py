from . import build_utils
from . import class_registry

@class_registry.register_class
class outputdefinition(class_registry.manage_state):
    '''
    outputdefinition Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.definitions = 'List of definitions'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Output definitions:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'definitions', 'List of potential outputs that can be requested, but which need additional data to be defined'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - outputdefinition Class'
        return s

