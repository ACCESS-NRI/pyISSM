import numpy as np
from . import param_utils
from . import class_registry
from .. import execute

## --------------------------------------------------------
## mesh.mesh2d
## --------------------------------------------------------
@class_registry.register_class
class mesh2d(class_registry.manage_state):
    """
    Horizontal 2D Mesh Class definition for ISSM.

    This class defines the default parameters for a 2D triangular mesh used in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    x : ndarray, default=np.nan
        Vertices x coordinate [m]
    y : ndarray, default=np.nan
        Vertices y coordinate [m]
    elements : ndarray, default=np.nan
        Vertex indices of the mesh elements
    numberofelements : int, default=0
        Number of elements
    numberofvertices : int, default=0
        Number of vertices
    numberofedges : int, default=0
        Number of edges of the 2d mesh
    lat : ndarray, default=np.nan
        Vertices latitude [degrees]
    long : ndarray, default=np.nan
        Vertices longitude [degrees]
    epsg : int, default=0
        EPSG code (ex: 3413 for UPS Greenland, 3031 for UPS Antarctica)
    scale_factor : float, default=np.nan
        Projection correction for volume, area, etc. computation
    vertexonboundary : ndarray, default=np.nan
        Vertices on the boundary of the domain flag list
    edges : ndarray, default=np.nan
        Edges of the 2d mesh (vertex1 vertex2 element1 element2)
    segments : ndarray, default=np.nan
        Edges on domain boundary (vertex1 vertex2 element)
    segmentmarkers : ndarray, default=np.nan
        Number associated to each segment
    vertexconnectivity : ndarray, default=np.nan
        List of elements connected to vertex_i
    elementconnectivity : ndarray, default=np.nan
        List of elements adjacent to element_i
    average_vertex_connectivity : int, default=25
        Average number of vertices connected to one vertex
    extractedvertices : ndarray, default=np.nan
        Vertices extracted from the model
    extractedelements : ndarray, default=np.nan
        Elements extracted from the model

    Methods
    -------
    __init__(self, other=None)
        Initializes the default mesh2d parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the mesh2d object.
    __str__(self)
        Returns a short string identifying the class.
    domain_type(self)
        Returns the domain type of the mesh.
    dimension(self)
        Returns the dimension of the mesh.
    element_type(self)
        Returns the element type of the mesh.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.mesh = pyissm.param.mesh.mesh2d()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.x = np.nan
        self.y = np.nan
        self.elements = np.nan
        self.numberofelements = 0
        self.numberofvertices = 0
        self.numberofedges = 0
        self.lat = np.nan
        self.long = np.nan
        self.epsg = 0
        self.scale_factor = np.nan
        self.vertexonboundary = np.nan
        self.edges = np.nan
        self.segments = np.nan
        self.segmentmarkers = np.nan
        self.vertexconnectivity = np.nan
        self.elementconnectivity = np.nan
        self.average_vertex_connectivity = 25
        self.extractedvertices = np.nan
        self.extractedelements = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   2D tria Mesh (horizontal):\n'

        s += '{}\n'.format('      Elements and vertices:')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'numberofelements', 'number of elements'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'numberofvertices', 'number of vertices'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'elements', 'vertex indices of the mesh elements'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'x', 'vertices x coordinate [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'y', 'vertices y coordinate [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'edges', 'edges of the 2d mesh (vertex1 vertex2 element1 element2)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'numberofedges', 'number of edges of the 2d mesh'))
        s += '\n'
        s += '{}\n'.format('      Properties:')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vertexonboundary', 'vertices on the boundary of the domain flag list'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'segments', 'edges on domain boundary (vertex1 vertex2 element)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'segmentmarkers', 'number associated to each segment'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vertexconnectivity', 'list of elements connected to vertex_i'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'elementconnectivity', 'list of elements adjacent to element_i'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'average_vertex_connectivity', 'average number of vertices connected to one vertex'))
        s += '\n'
        s += '{}\n'.format('      Extracted model:')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'extractedvertices', 'vertices extracted from the model'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'extractedelements', 'elements extracted from the model'))
        s += '\n'
        s += '{}\n'.format('      Projection:')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'lat', 'vertices latitude [degrees]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'long', 'vertices longitude [degrees]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'epsg', 'EPSG code (ex: 3413 for UPS Greenland, 3031 for UPS Antarctica)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'scale_factor', 'Projection correction for volume, area, etc. computation'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - mesh.mesh2d Class'
        return s

    # Define domain type
    def domain_type(self):
        return '2Dhorizontal'
    
    # Define mesh dimension
    def dimension(self):
        return 2

    # Define mesh element type
    def element_type(self):
        return 'Tria'

    # Marshall method for saving the mesh.mesh2d parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [mesh.mesh2d] parameters to a binary file.

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

        ## Write headers to file
        execute.WriteData(fid, prefix, name = 'md.mesh.domain_type', data = 'Domain' + self.domain_type(), format = 'String')
        execute.WriteData(fid, prefix, name = 'md.mesh.domain_dimension', data = self.dimension(), format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.mesh.elementtype', data = self.element_type(), format = 'String')

        ## Write Integer fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'numberofelements', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'numberofvertices', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'average_vertex_connectivity', format = 'Integer')

        ## Write DoubleMat fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'x', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'y', format = 'DoubleMat', mattype = 1 )
        execute.WriteData(fid, prefix, name = 'md.mesh.z', data = np.zeros(self.numberofvertices), format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'elements', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'vertexonboundary', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'segments', format = 'DoubleMat', mattype = 3)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'scale_factor', format = 'DoubleMat', mattype = 1)
        
        ## Write conditional fields
        if md.transient.isoceancoupling:
            execute.WriteData(fid, prefix, obj = self, fieldname = 'lat', format = 'DoubleMat', mattype = 1)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'long', format = 'DoubleMat', mattype = 1)

## --------------------------------------------------------
## mesh.mesh2dvertical
## --------------------------------------------------------
@class_registry.register_class
class mesh2dvertical(class_registry.manage_state):
    """
    Vertical 2D Mesh Class definition for ISSM.

    This class defines the default parameters for a 2D triangular mesh in the vertical plane used in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    x : ndarray, default=np.nan
        Vertices x coordinate [m]
    y : ndarray, default=np.nan
        Vertices y coordinate [m]
    elements : ndarray, default=np.nan
        Vertex indices of the mesh elements
    numberofelements : int, default=0
        Number of elements
    numberofvertices : int, default=0
        Number of vertices
    numberofedges : int, default=0
        Number of edges of the 2d mesh
    lat : ndarray, default=np.nan
        Vertices latitude [degrees]
    long : ndarray, default=np.nan
        Vertices longitude [degrees]
    epsg : int or float, default=np.nan
        EPSG code (ex: 3413 for UPS Greenland, 3031 for UPS Antarctica)
    scale_factor : float, default=np.nan
        Projection correction for volume, area, etc. computation
    vertexonboundary : ndarray, default=np.nan
        Vertices on the boundary of the domain flag list
    vertexonbase : ndarray, default=np.nan
        Vertices on the bed of the domain flag list
    vertexonsurface : ndarray, default=np.nan
        Vertices on the surface of the domain flag list
    edges : ndarray, default=np.nan
        Edges of the 2d mesh (vertex1 vertex2 element1 element2)
    segments : ndarray, default=np.nan
        Edges on domain boundary (vertex1 vertex2 element)
    segmentmarkers : ndarray, default=np.nan
        Number associated to each segment
    vertexconnectivity : ndarray, default=np.nan
        List of elements connected to vertex_i
    elementconnectivity : ndarray, default=np.nan
        List of elements adjacent to element_i
    average_vertex_connectivity : int, default=25
        Average number of vertices connected to one vertex

    Methods
    -------
    __init__(self, other=None)
        Initializes the default mesh2dvertical parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the mesh2dvertical object.
    __str__(self)
        Returns a short string identifying the class.
    domain_type(self)
        Returns the domain type of the mesh.
    dimension(self)
        Returns the dimension of the mesh.
    element_type(self)
        Returns the element type of the mesh.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.mesh = pyissm.param.mesh.mesh2dvertical()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.x = np.nan
        self.y = np.nan
        self.elements = np.nan
        self.numberofelements = 0
        self.numberofvertices = 0
        self.numberofedges = 0
        self.lat = np.nan
        self.long = np.nan
        self.epsg = np.nan
        self.scale_factor = np.nan
        self.vertexonboundary = np.nan
        self.vertexonbase = np.nan
        self.vertexonsurface = np.nan
        self.edges = np.nan
        self.segments = np.nan
        self.segmentmarkers = np.nan
        self.vertexconnectivity = np.nan
        self.elementconnectivity = np.nan
        self.average_vertex_connectivity = 25

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   2D tria Mesh (vertical):\n'

        s += '{}\n'.format('      Elements and vertices:')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'numberofelements', 'number of elements'))
        s += '{}\n'.format(param_utils.fielddisplay(self, "numberofvertices", "number of vertices"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "elements", "vertex indices of the mesh elements"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "x", "vertices x coordinate [m]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "y", "vertices y coordinate [m]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "edges", "edges of the 2d mesh (vertex1 vertex2 element1 element2)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "numberofedges", "number of edges of the 2d mesh"))
        s += '\n'
        s += '{}\n'.format('      Properties:')
        s += '{}\n'.format(param_utils.fielddisplay(self, "vertexonboundary", "vertices on the boundary of the domain flag list"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vertexonbase', 'vertices on the bed of the domain flag list'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vertexonsurface', 'vertices on the surface of the domain flag list'))
        s += '{}\n'.format(param_utils.fielddisplay(self, "segments", "edges on domain boundary (vertex1 vertex2 element)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "segmentmarkers", "number associated to each segment"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "vertexconnectivity", "list of elements connected to vertex_i"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "elementconnectivity", "list of elements adjacent to element_i"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "average_vertex_connectivity", "average number of vertices connected to one vertex"))
        s += '\n'
        s += '{}\n'.format('      Projection:')
        s += '{}\n'.format(param_utils.fielddisplay(self, "lat", "vertices latitude [degrees]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "long", "vertices longitude [degrees]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "epsg", "EPSG code (ex: 3413 for UPS Greenland, 3031 for UPS Antarctica)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "scale_factor", "Projection correction for volume, area, etc. computation"))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - mesh.mesh2dvertical Class'
        return s

    # Define domain type
    def domain_type(self):
        return '2Dvertical'
    
    # Define mesh dimension
    def dimension(self):
        return 2

    # Define mesh element type
    def element_type(self):
        return 'Tria'
    
    # Marshall method for saving the mesh.mesh2dvertical parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [mesh.mesh2dvertical] parameters to a binary file.

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

        ## Write headers to file
        execute.WriteData(fid, prefix, name = 'md.mesh.domain_type', data = 'Domain' + self.domain_type(), format = 'String')
        execute.WriteData(fid, prefix, name = 'md.mesh.domain_dimension', data = self.dimension(), format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.mesh.elementtype', data = self.element_type(), format = 'String')

        ## Write Integer fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'numberofelements', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'numberofvertices', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'average_vertex_connectivity', format = 'Integer')

        ## Write DoubleMat & BooleanMat fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'x', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'y', format = 'DoubleMat', mattype = 1 )
        execute.WriteData(fid, prefix, name = 'md.mesh.z', data = np.zeros(self.numberofvertices), format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'elements', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'vertexonbase', format = 'BooleanMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'vertexonsurface', format = 'BooleanMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'scale_factor', format = 'DoubleMat', mattype = 1)

## --------------------------------------------------------
## mesh.mesh3dprisms
## --------------------------------------------------------
@class_registry.register_class
class mesh3dprisms(class_registry.manage_state):
    """
    3D Prism Mesh Class definition for ISSM.

    This class defines the default parameters for a 3D mesh composed of extruded prisms, used in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    x : ndarray, default=np.nan
        Vertices x coordinate [m]
    y : ndarray, default=np.nan
        Vertices y coordinate [m]
    z : ndarray, default=np.nan
        Vertices z coordinate [m]
    elements : ndarray, default=np.nan
        Vertex indices of the 3D mesh elements
    numberoflayers : int, default=0
        Number of extrusion layers
    numberofelements : int, default=0
        Number of elements in the 3D mesh
    numberofvertices : int, default=0
        Number of vertices in the 3D mesh
    lat : ndarray, default=np.nan
        Vertices latitude [degrees]
    long : ndarray, default=np.nan
        Vertices longitude [degrees]
    epsg : int, default=0
        EPSG code (ex: 3413 for UPS Greenland, 3031 for UPS Antarctica)
    scale_factor : float, default=np.nan
        Projection correction for volume, area, etc. computation
    vertexonbase : ndarray, default=np.nan
        Lower vertices flags list
    vertexonsurface : ndarray, default=np.nan
        Upper vertices flags list
    lowerelements : ndarray, default=np.nan
        Lower element list (NaN for element on the lower layer)
    lowervertex : ndarray, default=np.nan
        Lower vertex list (NaN for vertex on the lower surface)
    upperelements : ndarray, default=np.nan
        Upper element list (NaN for element on the upper layer)
    uppervertex : ndarray, default=np.nan
        Upper vertex list (NaN for vertex on the upper surface)
    vertexonboundary : ndarray, default=np.nan
        Vertices on the boundary of the domain flag list
    vertexconnectivity : ndarray, default=np.nan
        List of elements connected to vertex_i
    elementconnectivity : ndarray, default=np.nan
        List of elements adjacent to element_i
    average_vertex_connectivity : int, default=25
        Average number of vertices connected to one vertex
    x2d : ndarray, default=np.nan
        2D mesh vertices x coordinate [m]
    y2d : ndarray, default=np.nan
        2D mesh vertices y coordinate [m]
    elements2d : ndarray, default=np.nan
        2D mesh vertex indices of the mesh elements
    numberofvertices2d : int, default=0
        Number of vertices in the original 2D mesh
    numberofelements2d : int, default=0
        Number of elements in the original 2D mesh
    extractedvertices : ndarray, default=np.nan
        Vertices extracted from the model
    extractedelements : ndarray, default=np.nan
        Elements extracted from the model

    Methods
    -------
    __init__(self, other=None)
        Initializes the default mesh3dprisms parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the mesh3dprisms object.
    __str__(self)
        Returns a short string identifying the class.
    domain_type(self)
        Returns the domain type of the mesh.
    dimension(self)
        Returns the dimension of the mesh.
    element_type(self)
        Returns the element type of the mesh.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.mesh = pyissm.param.mesh.mesh3dprisms()
    """

    # Initialise with default parameters
    def __init__(self, other = None):  # {{{
        self.x = np.nan
        self.y = np.nan
        self.z = np.nan
        self.elements = np.nan
        self.numberoflayers = 0
        self.numberofelements = 0
        self.numberofvertices = 0
        self.lat = np.nan
        self.long = np.nan
        self.epsg = 0
        self.scale_factor = np.nan
        self.vertexonbase = np.nan
        self.vertexonsurface = np.nan
        self.lowerelements = np.nan
        self.lowervertex = np.nan
        self.upperelements = np.nan
        self.uppervertex = np.nan
        self.vertexonboundary = np.nan
        self.vertexconnectivity = np.nan
        self.elementconnectivity = np.nan
        self.average_vertex_connectivity = 25
        self.x2d = np.nan
        self.y2d = np.nan
        self.elements2d = np.nan
        self.numberofvertices2d = 0
        self.numberofelements2d = 0
        self.extractedvertices = np.nan
        self.extractedelements = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):  # {{{
        s = '   3D prism Mesh:\n'

        s += '{}\n'.format('      Elements and vertices of the original 2d mesh3dprisms:')
        s += '{}\n'.format(param_utils.fielddisplay(self, "numberofelements2d", "number of elements"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "numberofvertices2d", "number of vertices"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "elements2d", "vertex indices of the mesh3dprisms elements"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "x2d", "vertices x coordinate [m]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "y2d", "vertices y coordinate [m]"))
        s += '\n'
        s += '{}\n'.format('      Elements and vertices of the extruded 3d mesh3dprisms:')
        s += '{}\n'.format(param_utils.fielddisplay(self, "numberofelements", "number of elements"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "numberofvertices", "number of vertices"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "elements", "vertex indices of the mesh3dprisms elements"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "x", "vertices x coordinate [m]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "y", "vertices y coordinate [m]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "z", "vertices z coordinate [m]"))
        s += '\n'
        s += '{}\n'.format('      Properties:')
        s += '{}\n'.format(param_utils.fielddisplay(self, "numberoflayers", "number of extrusion layers"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "vertexonbase", "lower vertices flags list"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "vertexonsurface", "upper vertices flags list"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "uppervertex", "upper vertex list (NaN for vertex on the upper surface)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "upperelements", "upper element list (NaN for element on the upper layer)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "lowervertex", "lower vertex list (NaN for vertex on the lower surface)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "lowerelements", "lower element list (NaN for element on the lower layer)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "vertexonboundary", "vertices on the boundary of the domain flag list"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "vertexconnectivity", "list of elements connected to vertex_i"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "elementconnectivity", "list of elements adjacent to element_i"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "average_vertex_connectivity", "average number of vertices connected to one vertex"))
        s += '\n'
        s += '{}\n'.format('      Extracted model:')
        s += '{}\n'.format(param_utils.fielddisplay(self, "extractedvertices", "vertices extracted from the model"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "extractedelements", "elements extracted from the model"))
        s += '\n'
        s += '{}\n'.format('      Projection:')
        s += '{}\n'.format(param_utils.fielddisplay(self, "lat", "vertices latitude [degrees]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "long", "vertices longitude [degrees]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "epsg", "EPSG code (ex: 3413 for UPS Greenland, 3031 for UPS Antarctica)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "scale_factor", "Projection correction for volume, area, etc. computation"))
        return s

        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - mesh.mesh3dprisms Class'
        return s

    # Define domain type
    def domain_type(self):
        return '3D'

    # Define mesh dimension
    def dimension(self):
        return 3

    # Define mesh element type
    def element_type(self):
        return 'Penta'

    # Marshall method for saving the mesh.mesh3dprisms parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [mesh.mesh3dprisms] parameters to a binary file.

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

        ## Write headers to file
        execute.WriteData(fid, prefix, name = 'md.mesh.domain_type', data = 'Domain' + self.domain_type(), format = 'String')
        execute.WriteData(fid, prefix, name = 'md.mesh.domain_dimension', data = self.dimension(), format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.mesh.elementtype', data = self.element_type(), format = 'String')

        ## Write Integer fields
        fieldnames = ['numberoflayers', 'numberofelements', 'numberofvertices', 'average_vertex_connectivity', 'numberofvertices2d', 'numberofelements2d']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Integer')

        ## Write DoubleMat & BooleanMat fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'x', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'y', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'z', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'elements', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'vertexonbase', format = 'BooleanMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'vertexonsurface', format = 'BooleanMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'lowerelements', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'upperelements', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'elements2d', format = 'DoubleMat', mattype = 3)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'scale_factor', format = 'DoubleMat', mattype = 1)

        ## Write conditional fields
        if md.transient.isoceancoupling:
            execute.WriteData(fid, prefix, obj = self, fieldname = 'lat', format = 'DoubleMat', mattype = 1)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'long', format = 'DoubleMat', mattype = 1)

## --------------------------------------------------------
## mesh.mesh3dsurface
## --------------------------------------------------------
@class_registry.register_class
class mesh3dsurface(class_registry.manage_state):
    """
    3D Surface Mesh Class definition for ISSM.

    This class defines the default parameters for a 3D triangular surface mesh used in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    x : ndarray, default=np.nan
        Vertices x coordinate [m]
    y : ndarray, default=np.nan
        Vertices y coordinate [m]
    z : ndarray, default=np.nan
        Vertices z coordinate [m]
    elements : ndarray, default=np.nan
        Vertex indices of the mesh elements
    numberofelements : int, default=0
        Number of elements
    numberofvertices : int, default=0
        Number of vertices
    numberofedges : int, default=0
        Number of edges of the 3d mesh
    lat : ndarray, default=np.nan
        Vertices latitude [degrees]
    long : ndarray, default=np.nan
        Vertices longitude [degrees]
    r : ndarray, default=np.nan
        Vertices radius [m]
    vertexonboundary : ndarray, default=np.nan
        Vertices on the boundary of the domain flag list
    edges : ndarray, default=np.nan
        Edges of the 3d mesh (vertex1 vertex2 element1 element2)
    segments : ndarray, default=np.nan
        Edges on domain boundary (vertex1 vertex2 element)
    segmentmarkers : ndarray, default=np.nan
        Number associated to each segment
    vertexconnectivity : ndarray, default=np.nan
        List of elements connected to vertex_i
    elementconnectivity : ndarray, default=np.nan
        List of elements adjacent to element_i
    average_vertex_connectivity : int, default=25
        Average number of vertices connected to one vertex
    extractedvertices : ndarray, default=np.nan
        Vertices extracted from the model
    extractedelements : ndarray, default=np.nan
        Elements extracted from the model

    Methods
    -------
    __init__(self, other=None)
        Initializes the default mesh3dsurface parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the mesh3dsurface object.
    __str__(self)
        Returns a short string identifying the class.
    domain_type(self)
        Returns the domain type of the mesh.
    dimension(self)
        Returns the dimension of the mesh.
    element_type(self)
        Returns the element type of the mesh.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.mesh = pyissm.param.mesh.mesh3dsurface()
    """

    # Initialise with default parameters
    def __init__(self, other = None):  # {{{
        self.x = np.nan
        self.y = np.nan
        self.z = np.nan
        self.elements = np.nan
        self.numberofelements = 0
        self.numberofvertices = 0
        self.numberofedges = 0
        self.lat = np.nan
        self.long = np.nan
        self.r = np.nan
        self.vertexonboundary = np.nan
        self.edges = np.nan
        self.segments = np.nan
        self.segmentmarkers = np.nan
        self.vertexconnectivity = np.nan
        self.elementconnectivity = np.nan
        self.average_vertex_connectivity = 25

        self.extractedvertices = np.nan
        self.extractedelements = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):  # {{{
        s = '   3D tria Mesh (surface):\n'

        s += '      Elements and vertices:'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'numberofelements', 'number of elements'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'numberofvertices', 'number of vertices'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'elements', 'vertex indices of the mesh elements'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'x', 'vertices x coordinate [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'y', 'vertices y coordinate [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'z', 'vertices z coordinate [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'lat', 'vertices latitude [degrees]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'long', 'vertices longitude [degrees]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'r', 'vertices radius [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'edges', 'edges of the 2d mesh (vertex1 vertex2 element1 element2)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'numberofedges', 'number of edges of the 2d mesh'))
        s +='\n'
        s += '      Properties:'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vertexonboundary', 'vertices on the boundary of the domain flag list'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'segments', 'edges on domain boundary (vertex1 vertex2 element)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'segmentmarkers', 'number associated to each segment'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vertexconnectivity', 'list of elements connected to vertex_i'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'elementconnectivity', 'list of elements adjacent to element_i'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'average_vertex_connectivity', 'average number of vertices connected to one vertex'))
        s += '\n'
        s += '      Extracted model():'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'extractedvertices', 'vertices extracted from the model()'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'extractedelements', 'elements extracted from the model()'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - mesh.mesh3dsurface Class'
        return s

    # Define domain type
    def domain_type(self):
        return '3Dsurface'
    
    # Define mesh dimension
    def dimension(self):
        return 2

    # Define mesh element type
    def element_type(self):
        return 'Tria'
    
    # Marshall method for saving the mesh.mesh3dsurface parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [mesh.mesh3dsurface] parameters to a binary file.

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

        ## Write headers to file
        execute.WriteData(fid, prefix, name = 'md.mesh.domain_type', data = 'Domain' + self.domain_type(), format = 'String')
        execute.WriteData(fid, prefix, name = 'md.mesh.domain_dimension', data = self.dimension(), format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.mesh.elementtype', data = self.element_type(), format = 'String')

        ## Write Integer fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'numberofelements', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'numberofvertices', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'average_vertex_connectivity', format = 'Integer')

        ## Write DoubleMat fields
        fieldnames = ['x', 'y', 'z', 'lat', 'long', 'r', 'vertexonboundary']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'DoubleMat', mattype = 1)
        
        execute.WriteData(fid, prefix, name = 'md.mesh.z', data = np.zeros(self.numberofvertices), format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'elements', format = 'DoubleMat', mattype = 2)