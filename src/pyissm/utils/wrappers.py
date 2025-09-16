"""
Python functions for C++ Python Wrappers in ISSM

This module contains various Python functions that point to Python wrappers for C++ modules in ISSM.

NOTE: Functionality here requires the following:
- ISSM has been installed using python wrappers.
- The Python version used to install ISSM Python wrappers is the same as that used to run pyISSM.
"""

import os
import importlib
import sys
import numpy as np
from .. import param

def load_issm_wrapper(func):
    """
    Decorator for ISSM Python wrapper functions that ensures the corresponding 
    compiled `_python` module exists and is loaded.

    This decorator:
    1. Checks that the environment variable `ISSM_DIR` is set.
    2. Adds the ISSM `lib` directory to `sys.path`.
    3. Checks that the shared library (`.so`) for the function exists.
    4. Imports the `_python` module.
    5. Attaches the `_python` function as an attribute (`_func`) to the wrapper.
    6. Returns the original wrapper function without calling it.

    Parameters
    ----------
    func : function
        The Python wrapper function to decorate. The decorator will attach 
        the corresponding compiled `_python` function as `func._func`.

    Returns
    -------
    function
        The original wrapper function with the `_python` function attached 
        as an attribute `_func`.

    Raises
    ------
    RuntimeError
        If `ISSM_DIR` is not set or if the corresponding `.so` file does not exist.

    Notes
    -----
    - The wrapper function itself is not executed when decorated.
    - The `_python` function can be called later via `func._func(*args, **kwargs)`.
    - This decorator assumes the `_python` module is named as `func.__name__ + '_python'`.
    """

    # 1. Ensure $ISSM_DIR is set
    issm_dir = os.environ.get("ISSM_DIR")
    if not issm_dir:
        raise RuntimeError(f"load_issm_wrapper: Environment variable ISSM_DIR is not set.\n\n"
                            "Ensure that ISSM is installed and the environment is properly configured.\n\n"
                            "add 'export ISSM_DIR=\"<path_to_issm_directory>\"'\n"
                            "     source $ISSM_DIR/etc/environment.sh\n\n"
                            "to your .bash_profile or .zprofile")

    # 1.1 Add to sys.path
    sys.path.append(os.path.join(issm_dir, "lib"))

    # 2. Determine module name and path
    module_name = f"{func.__name__}_python"
    so_path = os.path.join(issm_dir, "lib", f"{module_name}.so")

    # 3. Check that the .so exists
    if not os.path.exists(so_path):
        raise RuntimeError(f"load_issm_wrapper: Shared library '{so_path}' does not exist")

    # 4. Import the module
    module = importlib.import_module(module_name)

    # 5. Attach the loaded _python function to the wrapper
    func._func = getattr(module, module_name)

    # 6. Return the original wrapper function without calling it
    return func


## Triangle_python
@load_issm_wrapper
def Triangle(domain_outline_filename,
             rifts_filename = None,
             area = 1000000):
    """
    Generate a triangular mesh from a domain outline file.

    Wrapper function for $ISSM_DIR/lib/Triangle_python.

    Parameters
    ----------
    domain_outline_filename : str
        Path to an Argus domain outline file defining the mesh boundary.
    rifts_filename : str, optional, default=None
        Path to an Argus rifts file defining internal rifts within the domain.
    area : float, default=1000000
        Maximum area allowed for any element of the resulting mesh.


    Returns
    -------
    index : ndarray of int
        Array defining the triangulation connectivity (element indices).
    x : ndarray of float
        X coordinates of the mesh nodes.
    y : ndarray of float
        Y coordinates of the mesh nodes.
    segments : ndarray of int
        Array of exterior segments defining the domain outline.
    segmentmarkers : ndarray of int
        Array of flags marking each segment (for boundary conditions or identifiers).

    Examples
    --------
    >>> index, x, y, segments, segmentmarkers = Triangle("domain.exp", rifts=None, area=1000.0)
    """
    # Handle optional argument
    if rifts_filename is None:
        rifts_filename = ''

    # Call the loaded _python function
    return Triangle._func(domain_outline_filename, rifts_filename, area)

## BamgConvertMesh_python
@load_issm_wrapper
def BamgConvertMesh(index,
                    x,
                    y):
    """
    Convert a 2D mesh to BAMG geometry and BAMG mesh format.

    Wrapper function for $ISSM_DIR/lib/BamgConvertMesh_python.

    Parameters
    ----------
    index : array_like
        Indices of the mesh elements.
    x : array_like
        X coordinates of the mesh nodes.
    y : array_like
        Y coordinates of the mesh nodes.

    Returns
    -------
    bamggeom : BamgGeom
        The generated BAMG geometry object.
    bamgmesh : BamgMesh
        The generated BAMG mesh object.

    Examples
    --------
    >>> bamggeom, bamgmesh = BamgConvertMesh(md.mesh.elements, md.mesh.x, md.mesh.y)
    """

    # Call the loaded _python function
    return BamgConvertMesh._func(index, x, y)

## BamgMesher_python
@load_issm_wrapper
def BamgMesher(bamgmesh,
               bamggeom,
               bamgoptions):
    """
    Generate a BAMG mesh from a BAMG mesh, geometry and options.

    Wrapper function for $ISSM_DIR/lib/BamgMesher_python.

    Parameters
    ----------
    bamgmesh : BamgMesh
        The BAMG mesh object to refine.
    bamggeom : BamgGeom
        The BAMG geometry object defining the domain.
    bamgoptions : BamgOptions
        The BAMG options object containing meshing parameters.

    Returns
    -------
    bamggeom : BamgGeom
        The generated BAMG geometry object.
    bamgmesh : BamgMesh
        The generated BAMG mesh object.

    Returns
    -------
    bamggeom : BamgGeom
        The generated BAMG geometry object.
    bamgmesh : BamgMesh
        The generated BAMG mesh object.

    Examples
    --------
    >>> bamggeom, bamgmesh = BamgMesher(bamgmesh, bamggeom, bamgoptions)
    """

    # Call the loaded _python function
    return BamgMesher._func(bamgmesh, bamggeom, bamgoptions)

## BamgTriangulate_python
@load_issm_wrapper
def BamgTriangulate(x, y):
    """
    Delaunay triangulation of a set of points.

    Wrapper function for $ISSM_DIR/lib/BamgMesher_python.

    Parameters
    ----------
    x : array_like
        X coordinates of the points.
    y : array_like
        Y coordinates of the points.

    Returns
    -------
    indices : array_like
        Indices of the triangles formed by the Delaunay triangulation.

    Examples
    --------
    >>> indices = BamgTriangulate(x, y)
    """

    # Call the loaded _python function
    return BamgTriangulate._func(x, y)

## MeshProfileIntersection_python
@load_issm_wrapper
def MeshProfileIntersection(index,
                            x,
                            y,
                            filename):
    """
    Intersection of mesh with profile segments from an Argus .exp file.

    Wrapper function for $ISSM_DIR/lib/MeshProfileIntersection_python.

    Takes a .exp file (made of several profiles) and figures out its 
    intersection with a triangular mesh.

    Parameters
    ----------
    index : array_like
        Indices of the mesh elements defining the triangulation connectivity.
    x : array_like
        X coordinates of the mesh nodes.
    y : array_like
        Y coordinates of the mesh nodes.
    filename : str
        Path to an Argus-style .exp file containing the segments. Can contain
        groups of disconnected segments.

    Returns
    -------
    segments : ndarray
        Array where each row represents a segment intersection with the mesh.
        Each row contains [x1, y1, x2, y2, element_id] where (x1, y1) and
        (x2, y2) are the segment extremities and element_id is the ID of the
        element containing the segment. The number of rows equals the number
        of segments intersecting the mesh.

    Examples
    --------
    >>> segments = MeshProfileIntersection(md.mesh.elements, md.mesh.x, md.mesh.y, "profiles.exp")
    """
    
    # Call the loaded _python function
    return MeshProfileIntersection._func(index, x, y, filename)

## CountourToMesh_python
@load_issm_wrapper
def ContourToMesh(index,
                  x,
                  y,
                  contour_name,
                  interp_type,
                  edge_value):
    
    """
    Flag the elements or nodes inside a contour.

    Wrapper function for $ISSM_DIR/lib/ContourToMesh_python.

    Parameters
    ----------
    index : array_like
        Mesh triangulation connectivity (element indices).
    x : array_like
        X coordinates of the mesh nodes.
    y : array_like
        Y coordinates of the mesh nodes.
    contour_name : str
        Name of .exp file containing the contours.
    interp_type : str
        Type of interpolation. Must be one of 'element', 'node', or 'element and node'.
    edge_value : int
        Value (0, 1 or 2) associated to the nodes on the edges of the polygons.

    Returns
    -------
    in_nod : ndarray or tuple
        If `interp_type` is 'node': vector of flags (0 or 1) of size nel.
        If `interp_type` is 'element and node': tuple containing (in_nod, in_elem).
        Otherwise: not returned.
    in_elem : ndarray
        If `interp_type` is 'element': vector of flags (0 or 1) of size nel.
        If `interp_type` is 'element and node': returned as second element of tuple.
        Otherwise: not returned.

    Examples
    --------
    >>> in_nod = ContourToMesh(md.elements, md.x, md.y, 'Contour.exp', 'node', 1)
    >>> in_elements = ContourToMesh(md.elements, md.x, md.y, 'Contour.exp', 'element', 0)
    >>> in_nodes, in_elements = ContourToMesh(md.elements, md.x, md.y, 'Contour.exp', 'element and node', 0)
    """

    # Call the loaded _python function
    in_nod, in_elem = ContourToMesh._func(index, x, y, contour_name, interp_type, edge_value)

    # Conditionally return the appropriate values
    if interp_type == 'element':
        return in_elem
    elif interp_type == 'node':
        return in_nod
    elif interp_type == 'element and node':
        return in_nod, in_elem
    else:
        raise TypeError('interpolation type "{}" not supported yet'.format(interp_type))
    
## IssmConfig_python
@load_issm_wrapper
def IssmConfig(string):
    """
    Get ISSM configuration value for a specified parameter.

    Wrapper function for $ISSM_DIR/lib/IssmConfig_python.

    Parameters
    ----------
    string : str
        The configuration parameter name to retrieve.

    Returns
    -------
    value
        The configuration value associated with the specified parameter string.

    Examples
    --------
    >>> value = IssmConfig('parameter_name')
    >>> print(value)
    """

    # Call the loaded _python function
    return IssmConfig._func(string)

## ElementConnectivity_python
@load_issm_wrapper
def ElementConnectivity(elements, node_connectivity):
    """
    Build element connectivity using node connectivity and elements.

    Wrapper function for $ISSM_DIR/lib/ElementConnectivity_python.

    Parameters
    ----------
    elements : array_like
        Mesh triangulation connectivity (element indices).
    node_connectivity : array_like
        Node connectivity information.

    Returns
    -------
    element_connectivity : ndarray
        Element connectivity matrix.

    Examples
    --------
    >>> element_connectivity = ElementConnectivity(elements, node_connectivity)
    """

    # Call the loaded _python function
    ## NOTE: Value returned from wrapper function is a tuple, the first element of which being the result we actually want
    return ElementConnectivity._func(elements, node_connectivity)[0]

## NodeConnectivity_python
@load_issm_wrapper
def NodeConnectivity(elements, num_nodes):
    """
    Build node connectivity using element connectivity and nodes.

    Wrapper function for $ISSM_DIR/lib/NodeConnectivity_python.

    Parameters
    ----------
    elements : array_like
        Mesh triangulation connectivity (element indices).
    num_nodes : int
        Number of nodes in the mesh.

    Returns
    -------
    node_connectivity : ndarray
        Node connectivity matrix.

    Examples
    --------
    >>> node_connectivity = NodeConnectivity(elements, num_nodes)
    """
    # Call the loaded _python function
    ## NOTE: Value returned from wrapper function is a tuple, the first element of which being the result we actually want
    return NodeConnectivity._func(elements, num_nodes)[0]

## ContourToNodes_python
@load_issm_wrapper
def ContourToNodes(x,
                   y,
                   contourname,
                   edgevalue):
    """
    Flag vertices inside contour.

    Wrapper function for $ISSM_DIR/lib/ContourToNodes_python.

    Parameters
    ----------
    x : array_like
        X coordinates of the nodes.
    y : array_like
        Y coordinates of the nodes.
    contourname : str
        Name of .exp/.shp file containing the contours, or resulting structure 
        from call to expread/shpread.
    edgevalue : int
        Value (0, 1 or 2) defining the value associated to the nodes on the 
        edges of the polygons.

    Returns
    -------
    flags : ndarray
        Vector of flags (0 or 1) of size equal to the number of nodes.

    Examples
    --------
    >>> flags = ContourToNodes(x, y, 'contour.exp', 1)
    """

    # Call the loaded _python function
    return np.squeeze(ContourToNodes._func(x, y, contourname, edgevalue))

## InterpFromGridToMesh_python
@load_issm_wrapper
def InterpFromGridToMesh(x,
                         y,
                         data,
                         x_mesh,
                         y_mesh,
                         default_value):
    """
    Interpolation from a grid onto a list of points.

    Wrapper function for $ISSM_DIR/lib/InterpFromGridToMesh_python.

    Parameters
    ----------
    x : array_like
        X coordinates of the grid data (must be in increasing order).
    y : array_like
        Y coordinates of the grid data (must be in increasing order).
    data : array_like
        Matrix holding the data to be interpolated onto the mesh.
    x_mesh : array_like
        X coordinates of the points onto which we interpolate.
    y_mesh : array_like
        Y coordinates of the points onto which we interpolate.
    default_value : float
        Default value to use for points outside the grid.

    Returns
    -------
    data_mesh : ndarray
        Vector of mesh interpolated data.

    Examples
    --------
    >>> data_mesh = InterpFromGridToMesh(x_grid, y_grid, Vel, md.mesh.x, md.mesh.y, 0)
    """

    # Call the loaded _python function
    return np.squeeze(InterpFromGridToMesh._func(x, y, data, x_mesh, y_mesh, default_value))

## InterpFromMesh2d_python
@load_issm_wrapper
def InterpFromMesh2d(index,
                     x,
                     y,
                     data,
                     x_prime,
                     y_prime,
                     default_value = None,
                     contourname = None):
    """
    Interpolate data from a 2D mesh onto a set of points.

    Wrapper function for $ISSM_DIR/lib/InterpFromMesh2d_python.

    Parameters
    ----------
    index : array_like
        Index of the mesh where data is defined (triangulation connectivity).
    x : array_like
        X coordinates of the nodes where data is defined.
    y : array_like
        Y coordinates of the nodes where data is defined.
    data : array_like
        Vector holding the data to be interpolated onto the points.
    x_prime : array_like
        X coordinates of the mesh vertices onto which we interpolate.
    y_prime : array_like
        Y coordinates of the mesh vertices onto which we interpolate.
    default_value : float or array_like, optional
        Default value(s) to use for interpolation. Can be a scalar or vector 
        of size len(x_prime).
    contourname : str, optional
        Name of contour file. Linear interpolation will happen on all x_prime, 
        y_prime inside the contour, default value will be adopted on the rest 
        of the mesh.

    Returns
    -------
    data_prime : ndarray
        Vector of interpolated data at the target points.

    Examples
    --------
    >>> data_prime = InterpFromMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, data, x_new, y_new)
    >>> data_prime = InterpFromMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, data, x_new, y_new, default_value=0.0)
    >>> data_prime = InterpFromMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, data, x_new, y_new, default_value=0.0, contourname='contour.exp')
    """

    # Call the loaded _python function
    if default_value is None and contourname is None:
        data_prime = InterpFromMesh2d._func(index, x, y, data, x_prime, y_prime)
    elif not default_value is None and contourname is None:
        data_prime = InterpFromMesh2d._func(index, x, y, data, x_prime, y_prime, default_value)
    elif not default_value is None and not contourname is None:
        data_prime = InterpFromMesh2d._func(index, x, y, data, x_prime, y_prime, default_value, contourname)
    else:
        raise Exception('utils.wrappers.InterpFromMesh2d:: When defining a contourname, default_value must also be defined.')

    ## NOTE: Value returned from wrapper function is a tuple, the first element of which being the result we actually want
    return data_prime[0]

## InterpFromMeshToGrid_python
@load_issm_wrapper
def InterpFromMeshToGrid(index,
                         x,
                         y,
                         data,
                         xgrid,
                         ygrid,
                         default_value):
    """
    Interpolate data from a mesh onto a regular grid.

    Wrapper function for $ISSM_DIR/lib/InterpFromMeshToGrid_python.

    Parameters
    ----------
    index : array_like
        Index of the mesh where data is defined (triangulation connectivity).
    x : array_like
        X coordinates of the nodes where data is defined.
    y : array_like
        Y coordinates of the nodes where data is defined.
    data : array_like
        Vertex values of data to be interpolated onto the grid.
    xgrid : array_like
        X coordinates defining the grid.
    ygrid : array_like
        Y coordinates defining the grid.
    default_value : float
        Value assigned to points located outside the mesh.

    Returns
    -------
    grid_data : ndarray
        Interpolated data on the regular grid.

    Examples
    --------
    >>> grid_data = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y, data, xgrid, ygrid, default_value=0.0)
    """

    # Call the loaded _python function
    return np.squeeze(InterpFromMeshToGrid._func(index, x, y, data, xgrid, ygrid, default_value))

## InterpFromMeshToMesh2d_python
@load_issm_wrapper
def InterpFromMeshToMesh2d(index,
                           x,
                           y,
                           data,
                           x_interp,
                           y_interp,
                           default_value = None):
    """
    Interpolate from a 2D triangular mesh onto a list of points.

    Wrapper function for $ISSM_DIR/lib/InterpFromMeshToMesh2d_python.

    Parameters
    ----------
    index : array_like
        Index of the mesh where data is defined (triangulation connectivity).
    x : array_like
        X coordinates of the nodes where data is defined.
    y : array_like
        Y coordinates of the nodes where data is defined.
    data : array_like
        Matrix holding the data to be interpolated onto the mesh (one column per field).
    x_interp : array_like
        X coordinates of the points onto which we interpolate.
    y_interp : array_like
        Y coordinates of the points onto which we interpolate.
    default_value : float, optional
        Default value if point is outside of triangulation (instead of linear interpolation).

    Returns
    -------
    data_prime : ndarray
        Vector of interpolated data at the target points.

    Examples
    --------
    >>> interpolated_temp = InterpFromMeshToMesh2d(index, x, y, temperature, md.mesh.x, md.mesh.y)
    >>> interpolated_temp = InterpFromMeshToMesh2d(index, x, y, temperature, md.mesh.x, md.mesh.y, default_value = 253)
    """

    # Call the loaded _python function
    if default_value is None:
        data_prime = InterpFromMeshToMesh2d._func(index, x, y, data, x_interp, y_interp)
    elif not default_value is None:
        data_prime = InterpFromMeshToMesh2d._func(index, x, y, data, x_interp, y_interp, default_value)
    else:
        raise Exception('utils.wrappers.InterpFromMeshToMesh2d:: Something went wrong! Make sure you have provided all required arguments.')

    ## NOTE: Value returned from wrapper function is a tuple, the first element of which being the result we actually want
    return data_prime[0]

## InterpFromMeshToMesh3d_python
@load_issm_wrapper
def InterpFromMeshToMesh3d(index,
                           x,
                           y,
                           z,
                           data,
                           x_prime,
                           y_prime,
                           z_prime,
                           default_value):
    """
    Interpolate from a 3D hexahedron mesh onto a list of points.

    Wrapper function for $ISSM_DIR/lib/InterpFromMeshToMesh3d_python.

    Parameters
    ----------
    index : array_like
        Index of the mesh where data is defined (hexahedron connectivity).
    x : array_like
        X coordinates of the nodes where data is defined.
    y : array_like
        Y coordinates of the nodes where data is defined.
    z : array_like
        Z coordinates of the nodes where data is defined.
    data : array_like
        Matrix holding the data to be interpolated onto the mesh.
    x_prime : array_like
        X coordinates of the points onto which we interpolate.
    y_prime : array_like
        Y coordinates of the points onto which we interpolate.
    z_prime : array_like
        Z coordinates of the points onto which we interpolate.
    default_value : float, optional
        Default value if no data is found (holes).

    Returns
    -------
    data_prime : ndarray
        Vector of interpolated data at the target points.

    Examples
    --------
    >>> interpolated_temp = InterpFromMeshToMesh3d(index, x, y, z, temperature, md.mesh.x, md.mesh.y, md.mesh.z, 253)
    """

    # Call the loaded _python function
    return InterpFromMeshToMesh3d._func(index, x, y, z, data, x_prime, y_prime, z_prime, default_value)

## MeshPartition_python
@load_issm_wrapper
def MeshPartition(md,
                  n_partitions):
    """
    Partition mesh according to the number of areas, using Metis library.

    Wrapper function for $ISSM_DIR/lib/MeshPartition_python.

    Parameters
    ----------
    md : object
        ISSM model object containing the mesh to be partitioned.
    n_partitions : int
        Number of partitions to divide the mesh into.

    Returns
    -------
    element_partitioning : ndarray
        Vector of partitioning area numbers, for every element.
    node_partitioning : ndarray
        Vector of partitioning area numbers, for every node.

    Examples
    --------
    >>> element_partitioning, node_partitioning = MeshPartition(md, 4)
    """

    # Get mesh info from md.mesh
    n_vertices = md.mesh.numberofvertices
    elements = md.mesh.elements
    n_vertices_2d = 0
    n_layers = 1
    elements_2d = []

    # Conditional handling for different mesh types
    if isinstance(md.mesh, param.mesh.mesh3dprisms):
        element_type = md.mesh.element_type()
        n_vertices_2d = md.mesh.numberofvertices2d
        n_layers = md.mesh.numberoflayers
        elements_2d = md.mesh.elements2d
    elif isinstance(md.mesh, param.mesh.mesh2d):
        element_type = md.mesh.element_type()
    elif isinstance(md.mesh, param.mesh.mesh2dvertical):
        element_type = md.mesh.element_type()

    # Call the loaded _python function
    [element_partitioning, node_partitioning] = MeshPartition._func(n_vertices, elements, n_vertices_2d, elements_2d, n_layers, element_type, n_partitions)

    return [element_partitioning, node_partitioning]

## ProcessRifts_python
@load_issm_wrapper
def ProcessRifts(index,
                 x,
                 y,
                 segments,
                 segmentmarkers):
    """
    Split a mesh where a rift (or fault) is present.

    Wrapper function for $ISSM_DIR/lib/ProcessRifts_python.

    Parameters
    ----------
    index : array_like
        Initial triangulation connectivity (element indices).
    x : array_like
        X coordinates of the initial mesh nodes.
    y : array_like
        Y coordinates of the initial mesh nodes.
    segments : array_like
        Initial array of exterior segments defining the domain outline.
    segmentmarkers : array_like
        Initial array of flags marking each segment.

    Returns
    -------
    index_prime : ndarray
        Resulting triangulation connectivity after rift processing.
    x_prime : ndarray
        X coordinates of the resulting mesh nodes.
    y_prime : ndarray
        Y coordinates of the resulting mesh nodes.
    segments_prime : ndarray
        Resulting array of exterior segments.
    segmentmarkers_prime : ndarray
        Resulting array of flags marking each segment.
    rifts : ndarray
        Array containing the processed rift information.

    Examples
    --------
    >>> index_prime, x_prime, y_prime, segments_prime, segmentmarkers_prime, rifts = ProcessRifts(index, x, y, segments, segmentmarkers)
    """
    
    # Call the loaded _python function
    index_prime, x_prime, y_prime, segments_prime, segmentmarkers_prime, rifts = ProcessRifts._func(index, x, y, segments, segmentmarkers)

    return index_prime, x_prime, y_prime, segments_prime, segmentmarkers_prime, rifts

## ExpToLevelSet_python
@load_issm_wrapper
def ExpToLevelSet(x, y, contourname):
    """
    Determine levelset distance between a contour and a cloud of points.

    Wrapper function for $ISSM_DIR/lib/ExpToLevelSet_python.

    Parameters
    ----------
    x : array_like
        X coordinates of the cloud points.
    y : array_like
        Y coordinates of the cloud points.
    contourname : str
        Name of .exp file containing the contours.

    Returns
    -------
    distance : ndarray
        Distance vector representing a levelset where the 0 level is one 
        of the contour segments.

    Examples
    --------
    >>> distance = ExpToLevelSet(md.mesh.x, md.mesh.y, 'Contour.exp')
    """
    
    # Call the loaded _python function
    return ExpToLevelSet._func(x, y, contourname)
