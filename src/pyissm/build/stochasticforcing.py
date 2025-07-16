import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class stochasticforcing(class_registry.manage_state):
    '''
    stochasticforcing Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.isstochasticforcing = 0
        self.fields = 'List of fields'
        self.defaultdimension = 0
        self.default_id = np.nan
        self.covariance = np.nan
        self.timecovariance = np.nan
        self.stochastictimestep = 0
        self.randomflag = 1

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   stochasticforcing parameters:\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isstochasticforcing', 'is stochasticity activated?'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'fields', 'fields with stochasticity applied, ex: [\'SMBautoregression\'], or [\'SMBforcing\',\'DefaultCalving\']'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'defaultdimension', 'dimensionality of the noise terms (does not apply to fields with their specific dimension)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'default_id', 'id of each element for partitioning of the noise terms (does not apply to fields with their specific partition)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'covariance', 'covariance matrix for within- and between-fields covariance (units must be squared field units),multiple matrices can be concatenated along 3rd dimension to apply different covariances in time'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'timecovariance', 'starting dates at which covariances apply (only applicabe if multiple covariance matrices are prescribed)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'stochastictimestep', 'timestep at which new stochastic noise terms are generated (default: md.timestepping.time_step)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'randomflag', 'whether to apply real randomness (true) or pseudo-randomness with fixed seed (false)'))
        s += 'Available fields:\n'
        s += '   BasalforcingsDeepwaterMeltingRatearma\n'
        s += '   BasalforcingsSpatialDeepwaterMeltingRate\n'
        s += '   DefaultCalving\n'
        s += '   FloatingMeltRate\n'
        s += '   FrictionWaterPressure\n'
        s += '   FrictionCoulombWaterPressure\n'
        s += '   FrictionSchoofWaterPressure\n'
        s += '   FrontalForcingsRignotarma (thermal forcing)\n'
        s += '   FrontalForcingsSubglacialDischargearma\n'
        s += '   SMBarma\n'
        s += '   SMBforcing\n'
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - stochasticforcing Class'
        return s

