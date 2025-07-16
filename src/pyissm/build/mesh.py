import numpy as np
from . import build_utils
from . import class_registry

## --------------------------------------------------------
## mesh.mesh2d
## --------------------------------------------------------
@class_registry.register_class
class mesh2d(class_registry.manage_state):
    '''
    mesh.mesh2d Class definition
    '''

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
        s += '{}\n'.format(build_utils.fielddisplay(self, 'numberofelements', 'number of elements'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'numberofvertices', 'number of vertices'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'elements', 'vertex indices of the mesh elements'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'x', 'vertices x coordinate [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'y', 'vertices y coordinate [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'edges', 'edges of the 2d mesh (vertex1 vertex2 element1 element2)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'numberofedges', 'number of edges of the 2d mesh'))
        s += '\n'
        s += '{}\n'.format('      Properties:')
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vertexonboundary', 'vertices on the boundary of the domain flag list'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'segments', 'edges on domain boundary (vertex1 vertex2 element)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'segmentmarkers', 'number associated to each segment'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vertexconnectivity', 'list of elements connected to vertex_i'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'elementconnectivity', 'list of elements adjacent to element_i'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'average_vertex_connectivity', 'average number of vertices connected to one vertex'))
        s += '\n'
        s += '{}\n'.format('      Extracted model:')
        s += '{}\n'.format(build_utils.fielddisplay(self, 'extractedvertices', 'vertices extracted from the model'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'extractedelements', 'elements extracted from the model'))
        s += '\n'
        s += '{}\n'.format('      Projection:')
        s += '{}\n'.format(build_utils.fielddisplay(self, 'lat', 'vertices latitude [degrees]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'long', 'vertices longitude [degrees]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'epsg', 'EPSG code (ex: 3413 for UPS Greenland, 3031 for UPS Antarctica)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'scale_factor', 'Projection correction for volume, area, etc. computation'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - mesh.mesh2d Class'
        return s

## --------------------------------------------------------
## mesh.mesh2dvertical
## --------------------------------------------------------
@class_registry.register_class
class mesh2dvertical(class_registry.manage_state):
    '''
    mesh.mesh2dvertical Class definition
    '''

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
        s += '{}\n'.format(build_utils.fielddisplay(self, 'numberofelements', 'number of elements'))
        s += '{}\n'.format(build_utils.fielddisplay(self, "numberofvertices", "number of vertices"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "elements", "vertex indices of the mesh elements"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "x", "vertices x coordinate [m]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "y", "vertices y coordinate [m]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "edges", "edges of the 2d mesh (vertex1 vertex2 element1 element2)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "numberofedges", "number of edges of the 2d mesh"))
        s += '\n'
        s += '{}\n'.format('      Properties:')
        s += '{}\n'.format(build_utils.fielddisplay(self, "vertexonboundary", "vertices on the boundary of the domain flag list"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vertexonbase', 'vertices on the bed of the domain flag list'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vertexonsurface', 'vertices on the surface of the domain flag list'))
        s += '{}\n'.format(build_utils.fielddisplay(self, "segments", "edges on domain boundary (vertex1 vertex2 element)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "segmentmarkers", "number associated to each segment"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "vertexconnectivity", "list of elements connected to vertex_i"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "elementconnectivity", "list of elements adjacent to element_i"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "average_vertex_connectivity", "average number of vertices connected to one vertex"))
        s += '\n'
        s += '{}\n'.format('      Projection:')
        s += '{}\n'.format(build_utils.fielddisplay(self, "lat", "vertices latitude [degrees]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "long", "vertices longitude [degrees]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "epsg", "EPSG code (ex: 3413 for UPS Greenland, 3031 for UPS Antarctica)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "scale_factor", "Projection correction for volume, area, etc. computation"))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - mesh.mesh2dvertical Class'
        return s

## --------------------------------------------------------
## mesh.mesh3dprisms
## --------------------------------------------------------
@class_registry.register_class
class mesh3dprisms(class_registry.manage_state):
    '''
    mesh.mesh3dprisms Class definition
    '''

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
        s += '{}\n'.format(build_utils.fielddisplay(self, "numberofelements2d", "number of elements"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "numberofvertices2d", "number of vertices"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "elements2d", "vertex indices of the mesh3dprisms elements"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "x2d", "vertices x coordinate [m]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "y2d", "vertices y coordinate [m]"))
        s += '\n'
        s += '{}\n'.format('      Elements and vertices of the extruded 3d mesh3dprisms:')
        s += '{}\n'.format(build_utils.fielddisplay(self, "numberofelements", "number of elements"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "numberofvertices", "number of vertices"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "elements", "vertex indices of the mesh3dprisms elements"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "x", "vertices x coordinate [m]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "y", "vertices y coordinate [m]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "z", "vertices z coordinate [m]"))
        s += '\n'
        s += '{}\n'.format('      Properties:')
        s += '{}\n'.format(build_utils.fielddisplay(self, "numberoflayers", "number of extrusion layers"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "vertexonbase", "lower vertices flags list"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "vertexonsurface", "upper vertices flags list"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "uppervertex", "upper vertex list (NaN for vertex on the upper surface)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "upperelements", "upper element list (NaN for element on the upper layer)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "lowervertex", "lower vertex list (NaN for vertex on the lower surface)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "lowerelements", "lower element list (NaN for element on the lower layer)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "vertexonboundary", "vertices on the boundary of the domain flag list"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "vertexconnectivity", "list of elements connected to vertex_i"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "elementconnectivity", "list of elements adjacent to element_i"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "average_vertex_connectivity", "average number of vertices connected to one vertex"))
        s += '\n'
        s += '{}\n'.format('      Extracted model:')
        s += '{}\n'.format(build_utils.fielddisplay(self, "extractedvertices", "vertices extracted from the model"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "extractedelements", "elements extracted from the model"))
        s += '\n'
        s += '{}\n'.format('      Projection:')
        s += '{}\n'.format(build_utils.fielddisplay(self, "lat", "vertices latitude [degrees]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "long", "vertices longitude [degrees]"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "epsg", "EPSG code (ex: 3413 for UPS Greenland, 3031 for UPS Antarctica)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, "scale_factor", "Projection correction for volume, area, etc. computation"))
        return s

        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - mesh.mesh3dprisms Class'
        return s


## --------------------------------------------------------
## mesh.mesh3dsurface
## --------------------------------------------------------
@class_registry.register_class
class mesh3dsurface(class_registry.manage_state):
    '''
    mesh.mesh3dsurface Class definition
    '''

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
        s += '{}\n'.format(build_utils.fielddisplay(self, 'numberofelements', 'number of elements'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'numberofvertices', 'number of vertices'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'elements', 'vertex indices of the mesh elements'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'x', 'vertices x coordinate [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'y', 'vertices y coordinate [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'z', 'vertices z coordinate [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'lat', 'vertices latitude [degrees]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'long', 'vertices longitude [degrees]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'r', 'vertices radius [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'edges', 'edges of the 2d mesh (vertex1 vertex2 element1 element2)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'numberofedges', 'number of edges of the 2d mesh'))
        s +='\n'
        s += '      Properties:'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vertexonboundary', 'vertices on the boundary of the domain flag list'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'segments', 'edges on domain boundary (vertex1 vertex2 element)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'segmentmarkers', 'number associated to each segment'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'vertexconnectivity', 'list of elements connected to vertex_i'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'elementconnectivity', 'list of elements adjacent to element_i'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'average_vertex_connectivity', 'average number of vertices connected to one vertex'))
        s += '\n'
        s += '      Extracted model():'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'extractedvertices', 'vertices extracted from the model()'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'extractedelements', 'elements extracted from the model()'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - mesh.mesh3dsurface Class'
        return s