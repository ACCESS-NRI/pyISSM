import numpy as np
from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm.model import execute

@class_registry.register_class
class nodalvalue(class_registry.manage_state):
    """
    Nodal value parameters class for ISSM.

    This class encapsulates parameters for extracting nodal values from ISSM (Ice Sheet System Model) simulations.
    It allows users to specify particular nodes or vertices from which to extract field values during the simulation,
    providing a way to monitor specific locations over time.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    name : str, default=''
        Identifier for this nodalvalue response.
    definitionstring : str, default=''
        String that identifies this output definition uniquely, from 'Outputdefinition[1-10]'.
    model_string : str, default=''
        String for field that is being retrieved.
    node : float, default=nan
        Vertex index at which we retrieve the value.

    Methods
    -------
    __init__(self, other=None)
        Initializes the nodalvalue parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the nodalvalue parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.nodalvalue = pyissm.model.classes.nodalvalue()
    md.nodalvalue.name = 'velocity_at_glacier_terminus'
    md.nodalvalue.model_string = 'Vel'
    md.nodalvalue.node = 1245
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.name = ''
        self.definitionstring = ''
        self.model_string = ''
        self.node = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   nodalvalue parameters:\n'

        s += '{}\n'.format(class_utils._field_display(self, 'name', 'identifier for this nodalvalue response'))
        s += '{}\n'.format(class_utils._field_display(self, 'definitionstring', 'string that identifies this output definition uniquely, from \'Outputdefinition[1-10]\''))
        s += '{}\n'.format(class_utils._field_display(self, 'model_string', 'string for field that is being retrieved'))
        s += '{}\n'.format(class_utils._field_display(self, 'node', 'vertex index at which we retrieve the value'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - nodalvalue Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        if not isinstance(self.name, str):
            raise Exception("pyissm.model.classes.nodalvalue.check_consistency: 'name' field should be a string!")
        OutputdefinitionStringArray = []
        for i in range(100):
            OutputdefinitionStringArray.append('Outputdefinition{}'.format(i))
        class_utils._check_field(md, fieldname = 'self.definitionstring', field = self.definitionstring, values = OutputdefinitionStringArray)
        class_utils._check_field(md, fieldname = 'self.node', field = self.node, values = range(md.mesh.numberofvertices))
        
        return md

    # Marshall method for saving the nodalvalue parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [nodalvalue] parameters to a binary file.

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

        ## Write fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'name', format = 'String')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'definitionstring', format = 'String')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'model_string', format = 'String')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'node', format = 'Integer')
