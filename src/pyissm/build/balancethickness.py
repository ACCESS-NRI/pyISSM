import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class balancethickness(class_registry.manage_state):
    """
    Balance thickness solution parameters class for ISSM.

    This class encapsulates parameters for the balance thickness solution in the ISSM (Ice Sheet System Model) framework.
    It allows users to configure thickness constraints, thickening rates, and stabilization parameters for 
    solving the balance thickness equation.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    spcthickness : float, default=nan
        Thickness constraints (NaN means no constraint) [m].
    thickening_rate : float, default=nan
        Ice thickening rate used in the mass conservation (dh/dt) [m/yr].
    stabilization : int, default=1
        Stabilization method: 0=None, 1=SU, 2=SSA's artificial diffusivity, 3=DG.
    omega : float, default=nan
        Stabilization parameter.
    slopex : float, default=nan
        Surface slope in x-direction for stabilization.
    slopey : float, default=nan
        Surface slope in y-direction for stabilization.

    Methods
    -------
    __init__(self, other=None)
        Initializes the balance thickness parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the balance thickness parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.balancethickness = pyissm.build.balancethickness()
    md.balancethickness.spcthickness = 100.0
    md.balancethickness.thickening_rate = 0.1
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.spcthickness = np.nan
        self.thickening_rate = np.nan
        self.stabilization = 1
        self.omega = np.nan
        self.slopex = np.nan
        self.slopey = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   balance thickness solution parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'spcthickness', 'thickness constraints (NaN means no constraint) [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thickening_rate', 'ice thickening rate used in the mass conservation (dh / dt) [m / yr]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'stabilization', "0: None, 1: SU, 2: SSA's artificial diffusivity, 3:DG"))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - balancethickness Class'
        return s

