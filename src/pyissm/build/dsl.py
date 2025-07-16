import numpy as np
from . import build_utils
from . import class_registry

## ------------------------------------------------------
## dsl.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    '''
    dsl.default Class definition
    '''

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
        s += '{}\n'.format(build_utils.fielddisplay(self, 'global_average_thermosteric_sea_level', 'Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m).'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sea_surface_height_above_geoid', 'Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m).'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sea_water_pressure_at_sea_floor', 'Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!).'))
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
    '''
    dsl.mme Class definition
    '''

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
        s += '{}\n'.format(build_utils.fielddisplay(self, 'modelid', 'index into the multi-model ensemble, determines which field will be used.'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'global_average_thermosteric_sea_level', 'Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m) for each ensemble.'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sea_surface_height_above_geoid', 'Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m) for each ensemble.'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sea_water_pressure_at_sea_floor', 'Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!) for each ensemble.'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - dsl mme Class'
        return s