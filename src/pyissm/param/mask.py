import numpy as np
from . import param_utils
from . import class_registry
from .. import execute

@class_registry.register_class
class mask(class_registry.manage_state):
    """
    Mask parameters class for ISSM.

    This class encapsulates mask parameters for the ISSM (Ice Sheet System Model) framework.
    It defines level-set functions that determine the presence of ice and ocean, allowing for
    precise tracking of ice fronts, coastlines, and grounding lines using signed distance functions.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    ice_levelset : ndarray, default=nan
        Level-set function for ice: presence of ice if < 0, icefront position if = 0, no ice if > 0.
    ocean_levelset : ndarray, default=nan
        Level-set function for ocean: presence of ocean if < 0, coastline/grounding line if = 0, no ocean if > 0.

    Methods
    -------
    __init__(self, other=None)
        Initializes the mask parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the mask parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.mask = pyissm.param.mask()
    md.mask.ice_levelset = ice_levelset_field
    md.mask.ocean_levelset = ocean_levelset_field
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.ice_levelset = np.nan
        self.ocean_levelset = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   mask parameters:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ice_levelset', 'presence of ice if < 0, icefront position if = 0, no ice if > 0'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ocean_levelset', 'presence of ocean if < 0, coastline/grounding line if = 0, no ocean if > 0'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - mask Class'
        return s

    # Marshall method for saving the mask parameters
    def marshall_class(self, prefix, md, fid):
        """
        Marshall the mask parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.

        Returns
        -------
        None
        """

        ## Write each field to the file (all fields are of the same type/format)
        fieldnames = list(self.__dict__.keys())
        for fieldname in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = fieldname, format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)