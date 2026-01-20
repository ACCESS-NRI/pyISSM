import numpy as np
import os
from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm.model import execute
from pyissm import tools

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
    md.massfluxatgate = pyissm.model.classes.massfluxatgate()
    md.massfluxatgate.name = 'terminus_flux'
    md.massfluxatgate.profilename = 'terminus_profile.shp'
    md.massfluxatgate.definitionstring = 'Outputdefinition1'
    """

    # Initialise with default parameters
    def __init__(self, other = None, **kwargs):
        self.name = ''
        self.definitionstring = ''
        self.profilename = ''
        self.segments = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   massfluxatgate parameters:\n'
        
        s += '{}\n'.format(class_utils.fielddisplay(self, 'name', 'identifier for this massfluxatgate response'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'definitionstring', 'string that identifies this output definition uniquely, from Outputdefinition[1 - 100]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'profilename', 'name of file (shapefile or argus file) defining a profile (or gate)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - massfluxatgate Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        if not isinstance(self.name, str):
            raise RuntimeError("pyissm.model.classes.massfluxatgate.check_consistency: 'name' field should be a string.")

        if not isinstance(self.profilename, str):
            raise RuntimeError("pyissm.model.classes.massfluxatgate.check_consistency: 'profilename' field should be a string.")

        allowed = [f'Outputdefinition{i}' for i in range(1, 100)]

        class_utils.check_field(md, field = self.definitionstring, values = allowed)

        # Check the profilename points to a file!:
        if not os.path.isfile(self.profilename):
            raise FileNotFoundError(f"Profile file for gate not found: {self.profilename}")

        return md
    
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
        if tools.wrappers.check_wrappers_installed():
            self.segments = tools.wrappers.MeshProfileIntersection(index = md.mesh.elements,
                                                                x = md.mesh.x,
                                                                y = md.mesh.y,
                                                                filename = self.profilename)[0]
        else:
            ## If wrappers are not installed, raise error as segments are required to marshall class
            raise RuntimeError('massfluxatgate.marshall_class: Python wrappers not installed. Unable to compute segments for mass flux variable, required to marshall class.')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'name', format = 'String')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'definitionstring', format = 'String')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'segments', format = 'DoubleMat', mattype = 1)

