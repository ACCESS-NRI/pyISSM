import numpy as np
from . import param_utils
from . import class_registry

## ------------------------------------------------------
## dsl.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    """
    Default DSL parameters class for ISSM.

    This class encapsulates the default parameters for dynamic sea level (DSL) in the ISSM (Ice Sheet System Model) framework.
    It defines the main DSL-related parameters.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    global_average_thermosteric_sea_level : float or ndarray, default=np.nan
        Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m).
    sea_surface_height_above_geoid : float or ndarray, default=np.nan
        Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m).
    sea_water_pressure_at_sea_floor : float or ndarray, default=np.nan
        Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!).

    Methods
    -------
    __init__(self, other=None)
        Initializes the DSL parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the DSL parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.dsl = pyissm.build.dsl.default()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.global_average_thermosteric_sea_level = np.nan
        self.sea_surface_height_above_geoid = np.nan
        self.sea_water_pressure_at_sea_floor = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   dsl parameters:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'global_average_thermosteric_sea_level', 'Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m).'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sea_surface_height_above_geoid', 'Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m).'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sea_water_pressure_at_sea_floor', 'Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!).'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - dsl Class'
        return s

## ------------------------------------------------------
## dsl.mme
## ------------------------------------------------------
@class_registry.register_class
class mme(class_registry.manage_state):
    """
    Multi-Model Ensemble (MME) DSL parameters class for ISSM.

    This class encapsulates the parameters for dynamic sea level (DSL) for a multi-model ensemble (MME) of CMIP5 outputs in the ISSM (Ice Sheet System Model) framework.
    It defines the main DSL-related parameters for each ensemble member.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    modelid : int, default=0
        Index into the multi-model ensemble, determines which field will be used.
    global_average_thermosteric_sea_level : float or ndarray, default=np.nan
        Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m) for each ensemble member.
    sea_surface_height_above_geoid : float or ndarray, default=np.nan
        Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m) for each ensemble member.
    sea_water_pressure_at_sea_floor : float or ndarray, default=np.nan
        Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!) for each ensemble member.

    Methods
    -------
    __init__(self, other=None)
        Initializes the MME DSL parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the MME DSL parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.dsl = pyissm.build.dsl.mme()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.modelid = 0
        self.global_average_thermosteric_sea_level = np.nan
        self.sea_surface_height_above_geoid = np.nan
        self.sea_water_pressure_at_sea_floor = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   dsl mme parameters:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'modelid', 'index into the multi-model ensemble, determines which field will be used.'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'global_average_thermosteric_sea_level', 'Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m) for each ensemble.'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sea_surface_height_above_geoid', 'Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m) for each ensemble.'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sea_water_pressure_at_sea_floor', 'Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!) for each ensemble.'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - dsl mme Class'
        return s