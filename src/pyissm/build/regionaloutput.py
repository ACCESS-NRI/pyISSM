import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class regionaloutput(class_registry.manage_state):
    '''
    regionaloutput Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.name = ''
        self.model= ''
        self.definitionstring = ''
        self.outputnamestring = ''
        self.mask = np.nan
        self.maskexpstring = ''

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   regionaloutput parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'name', 'identifier for this regional response'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'definitionstring', 'string that identifies this output definition uniquely, from Outputdefinition[1 - 100]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'outputnamestring', 'string that identifies the type of output you want, eg. IceVolume, TotalSmb, GroudedArea'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'mask', 'mask vectorial field which identifies the region of interest (value > 0 will be included)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'maskexpstring', 'name of Argus file that can be passed in to define the regional mask'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - regionaloutput Class'
        return s

