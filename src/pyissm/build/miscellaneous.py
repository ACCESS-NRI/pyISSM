from . import build_utils
from . import class_registry

@class_registry.register_class
class miscellaneous(class_registry.manage_state):
    '''
    miscellaneous Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.notes = ''
        self.name = ''
        self.dummy = 'Placeholder for dummy fields'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   miscellaneous parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'notes', 'notes in a cell of strings'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'name', 'model name'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'dummy', 'empty field to store some data'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - miscellaneous Class'
        return s

