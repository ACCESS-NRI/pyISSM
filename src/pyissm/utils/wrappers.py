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
    >>> bamggeom, bamgmesh = BamgConvertMesh_python(md.mesh.elements, md.mesh.x, md.mesh.y)
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
    >>> bamggeom, bamgmesh = BamgMesher_python(bamgmesh, bamggeom, bamgoptions)
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

# CountourToMesh_python
@load_issm_wrapper
def ContourToMesh(index,
                  x,
                  y,
                  contour_name,
                  interp_type,
                  edge_value):
    
    """
    Flag the elements or nodes inside a contour.

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

    return IssmConfig._func(string)
