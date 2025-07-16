import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class misfit(class_registry.manage_state):
    '''
    misfit Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.name = ''
        self.definitionstring = ''
        self.model_string = ''
        self.observation = np.nan
        self.observation_string = ''
        self.timeinterpolation = 'nearestneighbor'
        self.local = 1
        self.weights = np.nan
        self.weights_string = ''
        self.cumulated = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Misfit:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'name', 'identifier for this misfit response'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'definitionstring', 'string that identifies this output definition uniquely, from "Outputdefinition[1 - 10]"'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'model_string', 'string for field that is modeled'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'observation', 'observed field that we compare the model against'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'observation_string', 'observation string'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'local', 'is the response local to the elements, or global? (default is 1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'timeinterpolation', 'interpolation routine used to interpolate misfit between two time steps (default is "nearestneighbor"'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'weights', 'weights (at vertices) to apply to the misfit'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'weights_string', 'string for weights for identification purposes'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - misfit Class'
        return s

