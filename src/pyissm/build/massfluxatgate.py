import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class massfluxatgate(class_registry.manage_state):
    '''
    massfluxatgate Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.name = ''
        self.definitionstring = ''
        self.profilename = ''
        self.segments = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Massfluxatgate:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'name', 'identifier for this massfluxatgate response'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'definitionstring', 'string that identifies this output definition uniquely, from Outputdefinition[1 - 100]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'profilename', 'name of file (shapefile or argus file) defining a profile (or gate)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - massfluxatgate Class'
        return s

