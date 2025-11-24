import numpy as np
from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry

@class_registry.register_class
class radaroverlay(class_registry.manage_state):
    """
    Radar overlay parameters class for ISSM.

    This class encapsulates parameters for radar overlay functionality in the ISSM (Ice Sheet System Model) framework.
    It stores radar power data and corresponding coordinates for visualization and analysis purposes,
    allowing integration of radar observations with ice sheet model results.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    pwr : ndarray, default=nan
        Radar power image (matrix).
    x : ndarray, default=nan
        Corresponding x coordinates [m].
    y : ndarray, default=nan
        Corresponding y coordinates [m].

    Methods
    -------
    __init__(self, other=None)
        Initializes the radaroverlay parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the radaroverlay parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.radaroverlay = pyissm.model.classes.radaroverlay()
    md.radaroverlay.pwr = radar_power_matrix
    md.radaroverlay.x = x_coordinates
    md.radaroverlay.y = y_coordinates
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.pwr = np.nan
        self.x = np.nan
        self.y = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   radaroverlay parameters:\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'pwr', 'radar power image (matrix)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'x', 'corresponding x coordinates [m]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'y', 'corresponding y coordinates [m]'))

        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - radaroverlay Class'
        return s

    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        # No checks
        return md
