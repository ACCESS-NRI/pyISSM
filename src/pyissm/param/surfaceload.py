import numpy as np
from . import param_utils
from . import class_registry

@class_registry.register_class
class surfaceload(class_registry.manage_state):
    """
    Surface load parameters class for ISSM.

    This class encapsulates parameters for surface loading in the ISSM (Ice Sheet System Model) framework.
    Surface loads include ice thickness changes, water height changes, and other surface loads (e.g., sediments)
    that affect solid Earth deformation and sea level calculations.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    icethicknesschange : ndarray, default=nan
        Thickness change: ice height equivalent [mIce/yr].
    waterheightchange : ndarray, default=nan
        Water height change: water height equivalent [mWater/yr].
    other : ndarray, default=nan
        Other loads (sediments) [kg/m^2/yr].

    Methods
    -------
    __init__(self, other=None)
        Initializes the surfaceload parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the surfaceload parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.surfaceload = pyissm.param.surfaceload()
    md.surfaceload.icethicknesschange = ice_thickness_change
    md.surfaceload.waterheightchange = water_height_change
    md.surfaceload.other = sediment_load
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.icethicknesschange = np.nan
        self.waterheightchange = np.nan
        self.other = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surfaceload:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'icethicknesschange', 'thickness change: ice height equivalent [mIce/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'waterheightchange', 'water height change: water height equivalent [mWater/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'other', 'other loads (sediments) [kg/m^2/yr]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - surfaceload Class'
        return s

