import numpy as np

from pyissm.execute import WriteData
from . import param_utils
from . import class_registry
from .. import execute
from .. import utils

@class_registry.register_class
class massfluxatgate(class_registry.manage_state):
    """
    Mass flux at gate parameters class for ISSM.

    This class encapsulates parameters for calculating mass flux through specified gates (profiles) in the ISSM (Ice Sheet System Model) framework.
    It allows users to define linear profiles or gates across which to calculate ice mass flux,
    providing a way to monitor ice discharge through specific cross-sections.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    name : str, default=''
        Identifier for this massfluxatgate response.
    definitionstring : str, default=''
        String that identifies this output definition uniquely, from Outputdefinition[1-100].
    profilename : str, default=''
        Name of file (shapefile or argus file) defining a profile (or gate).
    segments : float, default=nan
        Segments defining the gate geometry. Generated internally from the profile file.

    Methods
    -------
    __init__(self, other=None)
        Initializes the massfluxatgate parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the massfluxatgate parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.massfluxatgate = pyissm.param.massfluxatgate()
    md.massfluxatgate.name = 'terminus_flux'
    md.massfluxatgate.profilename = 'terminus_profile.shp'
    md.massfluxatgate.definitionstring = 'Outputdefinition1'
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.name = ''
        self.definitionstring = ''
        self.profilename = ''
        self.segments = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   massfluxatgate parameters:\n'
        
        s += '{}\n'.format(param_utils.fielddisplay(self, 'name', 'identifier for this massfluxatgate response'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'definitionstring', 'string that identifies this output definition uniquely, from Outputdefinition[1 - 100]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'profilename', 'name of file (shapefile or argus file) defining a profile (or gate)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - massfluxatgate Class'
        return s
    
    # Marshall method for saving the massfluxatgate parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [massfluxatgate] parameters to a binary file.

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

        ## Create segments from the profilename
        self.segments = utils.wrappers.MeshProfileIntersection(index = md.mesh.elements,
                                                               x = md.mesh.x,
                                                               y = md.mesh.y,
                                                               profile_name = self.profilename)[0]

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'name', format = 'String')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'definitionstring', format = 'String')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'segments', format = 'DoubleMat', mattype = 1)

