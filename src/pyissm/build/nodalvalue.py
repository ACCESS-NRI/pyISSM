import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class nodalvalue(class_registry.manage_state):
    '''
    nodalvalue Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.name = ''
        self.definitionstring = ''
        self.model_string = ''
        self.node = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   nodalvalue parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'name', 'identifier for this nodalvalue response'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'definitionstring', 'string that identifies this output definition uniquely, from \'Outputdefinition[1-10]\''))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'model_string', 'string for field that is being retrieved'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'node', 'vertex index at which we retrieve the value'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - nodalvalue Class'
        return s

