import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class initialization(class_registry.manage_state):
    """
    Initialization field values class for ISSM.

    This class encapsulates initial field values for various physical quantities in the ISSM (Ice Sheet System Model) framework.
    It provides storage for initial conditions of velocity, pressure, temperature, and other state variables
    that are used to initialize ice sheet model simulations.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    vx : ndarray, default=nan
        x component of velocity [m/yr].
    vy : ndarray, default=nan
        y component of velocity [m/yr].
    vz : ndarray, default=nan
        z component of velocity [m/yr].
    vel : ndarray, default=nan
        velocity norm [m/yr].
    pressure : ndarray, default=nan
        pressure [Pa].
    temperature : ndarray, default=nan
        temperature [K].
    enthalpy : ndarray, default=nan
        enthalpy [J].
    waterfraction : ndarray, default=nan
        fraction of water in the ice.
    sediment_head : ndarray, default=nan
        sediment water head of subglacial system [m].
    epl_head : ndarray, default=nan
        epl water head of subglacial system [m].
    epl_thickness : ndarray, default=nan
        thickness of the epl [m].
    watercolumn : ndarray, default=nan
        thickness of subglacial water [m].
    hydraulic_potential : ndarray, default=nan
        Hydraulic potential (for GlaDS) [Pa].
    channelarea : ndarray, default=nan
        subglacial water channel area (for GlaDS) [m2].
    sealevel : ndarray, default=nan
        sea level [m].
    bottompressure : ndarray, default=nan
        bottom pressure [Pa].
    dsl : ndarray, default=nan
        dynamic sea level [m].
    str : ndarray, default=nan
        surface temperature rate [K/yr].
    sample : ndarray, default=nan
        Realization of a Gaussian random field.
    debris : ndarray, default=nan
        Surface debris layer [m].
    age : ndarray, default=nan
        Initial age [yr].

    Methods
    -------
    __init__(self, other=None)
        Initializes the initialization parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the initialization parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.initialization = pyissm.build.initialization()
    md.initialization.vx = vx_initial
    md.initialization.vy = vy_initial
    md.initialization.temperature = temp_initial
    """

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

