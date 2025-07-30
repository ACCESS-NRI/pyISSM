from . import param_utils
from . import class_registry
import numpy as np
from .. import execute

@class_registry.register_class
class rifts(class_registry.manage_state):
    """
    Rifts parameters class for ISSM.

    This class encapsulates parameters for modeling rifts in the ISSM (Ice Sheet System Model) framework.
    Rifts are fractures or cracks in ice sheets that can affect ice dynamics and calving processes.
    This class stores structural information about rifts including their geometry, properties, and melange characteristics.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    riftstruct : str, default='Rift structure'
        Structure containing all rift information (vertices coordinates, segments, type of melange, etc.).
    riftproperties : str, default='Rift properties'
        Rift properties including physical and mechanical characteristics.

    Methods
    -------
    __init__(self, other=None)
        Initializes the rifts parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the rifts parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.rifts = pyissm.param.rifts()
    md.rifts.riftstruct = rift_structure_data
    md.rifts.riftproperties = rift_properties_data
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.riftstruct = 'Rift structure'
        self.riftproperties = 'Rift properties'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   rift parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'riftstruct', 'structure containing all rift information (vertices coordinates, segments, type of melange, ...)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'riftproperties', 'rift properties'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - rifts Class'
        return s

    # Marshall method for saving the rifts parameters
    def marshall_class(self, prefix, md, fid):
        """
        Marshall the rifts parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.

        Returns
        -------
        None
        """
        # TODO: Implement marshalling logic for riftstruct and riftproperties. Set to 0 for now to pass errors at runtime.
        execute.WriteData(fid, prefix, name = 'md.rifts.numrifts', data = 0, format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.rifts.riftstruct', data = np.zeros((0, 12)), format = 'DoubleMat', mattype = 3)