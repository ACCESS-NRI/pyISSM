from . import param_utils
from . import class_registry

@class_registry.register_class
class miscellaneous(class_registry.manage_state):
    """
    Miscellaneous parameters class for ISSM.

    This class encapsulates miscellaneous parameters and metadata for the ISSM (Ice Sheet System Model) framework.
    It provides storage for model notes, names, and other auxiliary information that doesn't fit into other categories.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    notes : str, default=''
        Notes in a cell of strings for documentation purposes.
    name : str, default=''
        Model name or identifier.
    dummy : str, default='Placeholder for dummy fields'
        Empty field to store some data or as placeholder.

    Methods
    -------
    __init__(self, other=None)
        Initializes the miscellaneous parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the miscellaneous parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.miscellaneous = pyissm.param.miscellaneous()
    md.miscellaneous.notes = 'Model run for Antarctic ice sheet'
    md.miscellaneous.name = 'Antarctica_2024'
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.notes = ''
        self.name = ''
        self.dummy = 'Placeholder for dummy fields'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   miscellaneous parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'notes', 'notes in a cell of strings'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'name', 'model name'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'dummy', 'empty field to store some data'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - miscellaneous Class'
        return s

