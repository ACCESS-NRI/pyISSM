import numpy as np
import os
from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm.model import execute, mesh
from pyissm import tools

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
    md.regionaloutput = pyissm.model.classes.regionaloutput()
    md.regionaloutput.name = 'west_antarctica_volume'
    md.regionaloutput.outputnamestring = 'IceVolume'
    md.regionaloutput.mask = west_antarctica_mask
    md.regionaloutput.definitionstring = 'Outputdefinition1'
    """

    # Initialise with default parameters
    def __init__(self, other = None, **kwargs):
        self.name = ''
        self.definitionstring = ''
        self.outputnamestring = ''
        self.mask = np.nan
        self.maskexpstring = None

        # Inherit matching fields from provided class
        super().__init__(other)

        # apply kwargs
        for key, value in kwargs.items():
            setattr(self, key, value)

    # Define repr
    def __repr__(self):
        s = '   regionaloutput parameters:\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'name', 'identifier for this regional response'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'definitionstring', 'string that identifies this output definition uniquely, from Outputdefinition[1 - 100]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'outputnamestring', 'string that identifies the type of output you want, eg. IceVolume, TotalSmb, GroudedArea'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'mask', 'mask vectorial field which identifies the region of interest (value > 0 will be included)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'maskexpstring', 'name of Argus file that can be passed in to define the regional mask'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - regionaloutput Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude regionaloutput fields to 3D mesh
        """
        self.mask = mesh.project_3d(md, vector = self.mask, type = 'node')

        return self

    # Check model consistency
    def check_consistency(self, md, solution, analyses):

        if not isinstance(self.name, str):
            raise RuntimeError("pyissm.model.classes.regionaloutput.check_consistency: 'name' field should be a string!")

        if not isinstance(self.outputnamestring, str):
            raise RuntimeError("pyissm.model.classes.regionaloutput.check_consistency: 'outputnamestring' field should be a string!")

        if self.maskexpstring is not None:
            if not os.path.isfile(self.maskexpstring):
                raise RuntimeError("pyissm.model.classes.regionaloutput.check_consistency: file name for mask exp does not point to a legitimate file on disk!")
            else:
                self.mask = tools.wrappers.ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y, self.maskexpstring, 'node', 1)

        OutputdefinitionStringArray = []
        for i in range(1, 100):
            x = 'Outputdefinition' + str(i)
            OutputdefinitionStringArray.append(x)

        class_utils.check_field(md, field = self.definitionstring, values = OutputdefinitionStringArray)
        class_utils.check_field(md, field = self.mask, size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        
        return md
    
    # Marshall method for saving the regionaloutput parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [regionaloutput] parameters to a binary file.

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

        ## Create mask from exp if maskexpstring is defined
        if self.maskexpstring is not None:
            if tools.wrappers.check_wrappers_installed():
                self.mask = tools.wrappers.ContourToMesh(index = md.mesh.elements,
                                                        x = md.mesh.x,
                                                        y = md.mesh.y,
                                                        contour_name = self.maskexpstring,
                                                        interp_type = 'node',
                                                        edge_value = 1)
            else:
                ## If wrappers are not installed, raise error as mask is required to marshall class
                raise RuntimeError('regionaloutput.marshall_class: Python wrappers not installed. Unable to compute regional mask, required to marshall class.')

        ## Write fields
        execute.WriteData(fid, prefix, data = self.name, name = 'md.regionaloutput.name', format = 'String')
        execute.WriteData(fid, prefix, data = self.definitionstring, name = 'md.regionaloutput.definitionstring', format = 'String')
        execute.WriteData(fid, prefix, data = self.outputnamestring, name = 'md.regionaloutput.outputnamestring', format = 'String')
        execute.WriteData(fid, prefix, data = self.mask, name = 'md.regionaloutput.mask', format = 'DoubleMat', mattype = 1)
