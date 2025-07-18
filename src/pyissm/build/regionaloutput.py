import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class regionaloutput(class_registry.manage_state):
    """
    Regional output parameters class for ISSM.

    This class encapsulates parameters for defining regional outputs in the ISSM (Ice Sheet System Model) framework.
    It allows users to extract integrated quantities (like ice volume, mass balance, grounded area) 
    over specific regions of interest defined by masks, providing regional analysis capabilities.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    name : str, default=''
        Identifier for this regional response.
    definitionstring : str, default=''
        String that identifies this output definition uniquely, from Outputdefinition[1-100].
    outputnamestring : str, default=''
        String that identifies the type of output you want, e.g. IceVolume, TotalSmb, GroundedArea.
    mask : ndarray, default=nan
        Mask vectorial field which identifies the region of interest (value > 0 will be included).
    maskexpstring : str, default=''
        Name of Argus file that can be passed in to define the regional mask.

    Methods
    -------
    __init__(self, other=None)
        Initializes the regionaloutput parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the regionaloutput parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.regionaloutput = pyissm.build.regionaloutput()
    md.regionaloutput.name = 'west_antarctica_volume'
    md.regionaloutput.outputnamestring = 'IceVolume'
    md.regionaloutput.mask = west_antarctica_mask
    md.regionaloutput.definitionstring = 'Outputdefinition1'
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.name = ''
        self.definitionstring = ''
        self.outputnamestring = ''
        self.mask = np.nan
        self.maskexpstring = ''

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   regionaloutput parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'name', 'identifier for this regional response'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'definitionstring', 'string that identifies this output definition uniquely, from Outputdefinition[1 - 100]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'outputnamestring', 'string that identifies the type of output you want, eg. IceVolume, TotalSmb, GroudedArea'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'mask', 'mask vectorial field which identifies the region of interest (value > 0 will be included)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'maskexpstring', 'name of Argus file that can be passed in to define the regional mask'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - regionaloutput Class'
        return s

