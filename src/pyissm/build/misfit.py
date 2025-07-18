import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class misfit(class_registry.manage_state):
    """
    Misfit parameters class for ISSM.

    This class encapsulates parameters for misfit calculations in the ISSM (Ice Sheet System Model) framework.
    Misfit functions measure the difference between model predictions and observations,
    and are essential for model validation, calibration, and inverse problem solutions.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    name : str, default=''
        Identifier for this misfit response.
    definitionstring : str, default=''
        String that identifies this output definition uniquely, from "Outputdefinition[1-10]".
    model_string : str, default=''
        String for field that is modeled.
    observation : ndarray, default=nan
        Observed field that we compare the model against.
    observation_string : str, default=''
        Observation string for identification purposes.
    timeinterpolation : str, default='nearestneighbor'
        Interpolation routine used to interpolate misfit between two time steps.
    local : int, default=1
        Is the response local to the elements, or global?
    weights : ndarray, default=nan
        Weights (at vertices) to apply to the misfit.
    weights_string : str, default=''
        String for weights for identification purposes.
    cumulated : float, default=nan
        Cumulated misfit value.

    Methods
    -------
    __init__(self, other=None)
        Initializes the misfit parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the misfit parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.misfit = pyissm.build.misfit()
    md.misfit.name = 'velocity_misfit'
    md.misfit.model_string = 'Vel'
    md.misfit.observation = observed_velocity
    md.misfit.weights = velocity_weights
    md.misfit.timeinterpolation = 'linear'
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.name = ''
        self.definitionstring = ''
        self.model_string = ''
        self.observation = np.nan
        self.observation_string = ''
        self.timeinterpolation = 'nearestneighbor'
        self.local = 1
        self.weights = np.nan
        self.weights_string = ''
        self.cumulated = None

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Misfit:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'name', 'identifier for this misfit response'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'definitionstring', 'string that identifies this output definition uniquely, from "Outputdefinition[1 - 10]"'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'model_string', 'string for field that is modeled'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'observation', 'observed field that we compare the model against'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'observation_string', 'observation string'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'local', 'is the response local to the elements, or global? (default is 1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'timeinterpolation', 'interpolation routine used to interpolate misfit between two time steps (default is "nearestneighbor"'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'weights', 'weights (at vertices) to apply to the misfit'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'weights_string', 'string for weights for identification purposes'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - misfit Class'
        return s

