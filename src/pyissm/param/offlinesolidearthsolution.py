import numpy as np
from . import param_utils
from . import class_registry

@class_registry.register_class
class offlinesolidearthsolution(class_registry.manage_state):
    """
    OfflineSolidEarthSolution class for ISSM.

    This class defines the default parameters for the offline solid-Earth solution used in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    displacementeast : float or ndarray, default=np.nan
        Solid-Earth Eastwards bedrock displacement time series (m)
    displacementnorth : float or ndarray, default=np.nan
        Solid-Earth Northwards bedrock displacement time series (m)
    displacementup : float or ndarray, default=np.nan
        Solid-Earth bedrock uplift time series (m)
    geoid : float or ndarray, default=np.nan
        Solid-Earth geoid time series (m)

    Methods
    -------
    __init__(self, other=None)
        Initializes the default parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the object.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.solidearth = pyissm.build.offlinesolidearthsolution()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.displacementeast = np.nan
        self.displacementnorth = np.nan
        self.displacementup = np.nan
        self.geoid = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    def __repr__(self):
        s = '         units for time series is (yr)\n       external: offlinesolidearth solution\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'displacementeast', 'solid-Earth Eastwards bedrock displacement time series (m)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'displacementnorth', 'solid-Earth Northwards bedrock displacement time series (m)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'displacementup', 'solid-Earth bedrock uplift time series (m)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'geoid', 'solid-Earth geoid time series (m)'))
        return s

    ## Define class string
    def __str__(self):
        s = 'ISSM - offlinesolidearthsolution Class'
        return s

