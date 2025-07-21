from . import param_utils
from . import class_registry

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
    md.outputdefinition = pyissm.build.outputdefinition()
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

