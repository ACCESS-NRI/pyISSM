import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class geometry(class_registry.manage_state):
    """
    Geometry parameters class for ISSM.

    This class encapsulates geometric parameters that define the ice sheet geometry in the ISSM (Ice Sheet System Model) framework.
    It stores elevation data for ice surface, thickness, base, and bed that are fundamental to ice sheet modeling.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    surface : ndarray, default=nan
        Ice upper surface elevation [m].
    thickness : ndarray, default=nan
        Ice thickness [m].
    base : ndarray, default=nan
        Ice base elevation [m].
    bed : ndarray, default=nan
        Bed elevation [m].
    hydrostatic_ratio : float, default=nan
        Hydrostatic ratio for floating ice.

    Methods
    -------
    __init__(self, other=None)
        Initializes the geometry parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the geometry parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.geometry = pyissm.build.geometry()
    md.geometry.surface = surface_elevation
    md.geometry.thickness = ice_thickness
    md.geometry.bed = bed_elevation
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.surface = np.nan
        self.thickness = np.nan
        self.base = np.nan
        self.bed = np.nan
        self.hydrostatic_ratio = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   geometry parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'surface', 'ice upper surface elevation [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thickness', 'ice thickness [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'base', 'ice base elevation [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'bed', 'bed elevation [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'hydrostatic_ratio', 'hydrostatic ratio for floating ice'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - geometry Class'
        return s

