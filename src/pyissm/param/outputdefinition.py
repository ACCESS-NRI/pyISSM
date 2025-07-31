from . import param_utils
from . import class_registry
from .. import execute

@class_registry.register_class
class outputdefinition(class_registry.manage_state):
    """
    Output definition parameters class for ISSM.

    This class encapsulates parameters for defining custom outputs in the ISSM (Ice Sheet System Model) framework.
    It allows users to specify additional outputs that can be requested during simulations,
    providing flexibility in extracting specific quantities or derived fields from the model.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    definitions : str, default='List of definitions'
        List of potential outputs that can be requested, but which need additional data to be defined.

    Methods
    -------
    __init__(self, other=None)
        Initializes the outputdefinition parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the outputdefinition parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.outputdefinition = pyissm.param.outputdefinition()
    md.outputdefinition.definitions = ['IceVolume', 'IceVolumeAboveFloatation', 'CustomOutput1']
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.definitions = 'List of definitions'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Output definitions:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'definitions', 'List of potential outputs that can be requested, but which need additional data to be defined'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - outputdefinition Class'
        return s

    # Marshall method for saving the outputdefinition parameters
    def marshall_class(self, prefix, md, fid):
        """
        Marshall the outputdefinition parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.

        Returns
        -------
        None
        """
        ## TODO: Iteratre through definitions and write them individually
        ## For now, just write the definitions list as a string array
        data = []
        execute.WriteData(fid, prefix, name = 'md.outputdefinition.list', data = data, format = 'StringArray')