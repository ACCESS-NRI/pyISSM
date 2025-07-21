import numpy as np
from . import param_utils
from . import class_registry

@class_registry.register_class
class esa(class_registry.manage_state):
    """
    Elastic solid earth adjustment (ESA) parameters class for ISSM.

    This class encapsulates parameters for elastic solid earth adjustment in the ISSM (Ice Sheet System Model) framework.
    ESA models the instantaneous elastic response of the solid Earth to ice loading changes,
    which is important for studying present-day ice sheet mass balance and sea level change.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    deltathickness : float, default=nan
        Thickness change: ice height equivalent [m].
    love_h : float, default=0.0
        Load Love number for radial displacement.
    love_l : float, default=0.0
        Load Love number for horizontal displacements.
    hemisphere : float, default=0.0
        North-south, East-west components of 2-D horizontal displacement vector: -1 south, 1 north.
    degacc : float, default=0.01
        Accuracy (default 0.01 deg) for numerical discretization of the Green's functions.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested (default: EsaUmotion).
    transitions : str, default='List of transitions'
        Indices into parts of the mesh that will be icecaps.

    Methods
    -------
    __init__(self, other=None)
        Initializes the esa parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the esa parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.esa = pyissm.build.esa()
    md.esa.deltathickness = thickness_change
    md.esa.love_h = 0.6
    md.esa.love_l = 0.1
    md.esa.degacc = 0.005
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.deltathickness = np.nan
        self.love_h = 0.
        self.love_l = 0.
        self.hemisphere = 0.
        self.degacc = 0.01
        self.requested_outputs = 'List of requested outputs'
        self.transitions = 'List of transitions'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   esa parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'deltathickness', 'thickness change: ice height equivalent [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'love_h', 'load Love number for radial displacement'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'love_l', 'load Love number for horizontal displaements'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'hemisphere', 'North-south, East-west components of 2-D horiz displacement vector:-1 south, 1 north'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'degacc', 'accuracy (default .01 deg) for numerical discretization of the Green''s functions'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'transitions', 'indices into parts of the mesh that will be icecaps'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested (default: EsaUmotion)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - esa Class'
        return s

