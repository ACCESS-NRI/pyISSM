import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class initialization(class_registry.manage_state):
    '''
    initialization Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.vx = np.nan
        self.vy = np.nan
        self.vz = np.nan
        self.vel = np.nan
        self.pressure = np.nan
        self.temperature = np.nan
        self.enthalpy = np.nan
        self.waterfraction = np.nan
        self.sediment_head = np.nan
        self.epl_head = np.nan
        self.epl_thickness = np.nan
        self.watercolumn = np.nan
        self.hydraulic_potential = np.nan
        self.channelarea = np.nan
        self.sealevel = np.nan
        self.bottompressure = np.nan
        self.dsl = np.nan
        self.str = np.nan
        self.sample = np.nan
        self.debris = np.nan
        self.age = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   initial field values:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'vx', 'x component of velocity [m/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vy', 'y component of velocity [m/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vz', 'z component of velocity [m/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vel', 'velocity norm [m/yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'pressure', 'pressure [Pa]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'temperature', 'temperature [K]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'enthalpy', 'enthalpy [J]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'waterfraction', 'fraction of water in the ice'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'watercolumn', 'thickness of subglacial water [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sediment_head', 'sediment water head of subglacial system [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'epl_head', 'epl water head of subglacial system [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'epl_thickness', 'thickness of the epl [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'hydraulic_potential', 'Hydraulic potential (for GlaDS) [Pa]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'channelarea', 'subglaciale water channel area (for GlaDS) [m2]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sample', 'Realization of a Gaussian random field'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'debris', 'Surface debris layer [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'age', 'Initial age [yr]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - initialization Class'
        return s

