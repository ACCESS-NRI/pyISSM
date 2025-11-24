import numpy as np
from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm.model import execute

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
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.misfit = pyissm.model.classes.misfit()
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

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Misfit:\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'name', 'identifier for this misfit response'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'definitionstring', 'string that identifies this output definition uniquely, from "Outputdefinition[1 - 10]"'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'model_string', 'string for field that is modeled'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'observation', 'observed field that we compare the model against'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'observation_string', 'observation string'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'local', 'is the response local to the elements, or global? (default is 1)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'timeinterpolation', 'interpolation routine used to interpolate misfit between two time steps (default is "nearestneighbor"'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'weights', 'weights (at vertices) to apply to the misfit'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'weights_string', 'string for weights for identification purposes'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - misfit Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        if type(self.name) != str:
            raise TypeError('pyissm.model.classes.misfit.check_consistency: "name" field should be a string!')

        OutputdefinitionStringArray = []
        for i in range(100):
            OutputdefinitionStringArray.append('Outputdefinition' + str(i))

        class_utils.check_field(md, fieldname = 'self.definitionstring', field = self.definitionstring, values = OutputdefinitionStringArray)
        if type(self.timeinterpolation) != str:
            raise TypeError('pyissm.model.classes.misfit.check_consistency: "timeinterpolation" field should be a string!')

        class_utils.check_field(md, fieldname = 'self.observation', field = self.observation, timeseries = True, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'self.timeinterpolation', field = self.timeinterpolation, values = ['nearestneighbor'])
        class_utils.check_field(md, fieldname = 'self.weights', field = self.weights, timeseries = True, allow_nan = False, allow_inf = False)

        return md

    # Marshall method for saving the misfit parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [misfit] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
            ISSM model object needed in some cases.

        Returns
        -------
        None
        """

        ## Write String fields
        fieldnames = ['name', 'definitionstring', 'model_string', 'observation_string', 'timeinterpolation',
                      'weights_string']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj=self, fieldname=field, format='String')

        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'observation', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'local', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'weights', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        

