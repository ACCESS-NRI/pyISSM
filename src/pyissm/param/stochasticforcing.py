import numpy as np
from . import param_utils
from . import class_registry

@class_registry.register_class
class stochasticforcing(class_registry.manage_state):
    """
    Stochastic forcing parameters class for ISSM.

    This class encapsulates parameters for stochastic forcing in the ISSM (Ice Sheet System Model) framework.
    It allows users to apply random forcing to various physical processes such as surface mass balance,
    basal melting, and calving, enabling uncertainty quantification and probabilistic modeling.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    isstochasticforcing : int, default=0
        Is stochasticity activated?
    fields : str, default='List of fields'
        Fields with stochasticity applied, e.g. ['SMBautoregression'], or ['SMBforcing','DefaultCalving'].
    defaultdimension : int, default=0
        Dimensionality of the noise terms (does not apply to fields with their specific dimension).
    default_id : ndarray, default=nan
        ID of each element for partitioning of the noise terms (does not apply to fields with their specific partition).
    covariance : float, default=nan
        Covariance matrix for within- and between-fields covariance (units must be squared field units), multiple matrices can be concatenated along 3rd dimension to apply different covariances in time.
    timecovariance : float, default=nan
        Starting dates at which covariances apply (only applicable if multiple covariance matrices are prescribed).
    stochastictimestep : float, default=0
        Timestep at which new stochastic noise terms are generated (default: md.timestepping.time_step).
    randomflag : int, default=1
        Whether to apply real randomness (true) or pseudo-randomness with fixed seed (false).

    Methods
    -------
    __init__(self, other=None)
        Initializes the stochasticforcing parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the stochasticforcing parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.stochasticforcing = pyissm.param.stochasticforcing()
    md.stochasticforcing.isstochasticforcing = 1
    md.stochasticforcing.fields = ['SMBforcing']
    md.stochasticforcing.defaultdimension = 1
    md.stochasticforcing.covariance = covariance_matrix
    """

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
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isstochasticforcing', 'is stochasticity activated?'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'fields', 'fields with stochasticity applied, ex: [\'SMBautoregression\'], or [\'SMBforcing\',\'DefaultCalving\']'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'defaultdimension', 'dimensionality of the noise terms (does not apply to fields with their specific dimension)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'default_id', 'id of each element for partitioning of the noise terms (does not apply to fields with their specific partition)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'covariance', 'covariance matrix for within- and between-fields covariance (units must be squared field units),multiple matrices can be concatenated along 3rd dimension to apply different covariances in time'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'timecovariance', 'starting dates at which covariances apply (only applicabe if multiple covariance matrices are prescribed)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'stochastictimestep', 'timestep at which new stochastic noise terms are generated (default: md.timestepping.time_step)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'randomflag', 'whether to apply real randomness (true) or pseudo-randomness with fixed seed (false)'))
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

