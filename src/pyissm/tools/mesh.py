"""
Tools for creating ISSM model meshes.
"""

import numpy as np
import warnings
import os
import collections
from .. import utils
from . import exp

def triangle(md,
             domain_name,
             resolution,
             rift_name = None):
    
    """
    Create a triangular mesh for an ISSM model using Triangle mesh generator.

    This function generates a triangular mesh based on a domain outline and optional
    rift constraints. It uses the Triangle mesh generator to create high-quality
    Delaunay triangulations with specified resolution constraints.

    Parameters
    ----------
    md : object
        ISSM model object containing the mesh structure to be populated.
    domain_name : str
        Path to the file containing the domain outline geometry.
    resolution : float
        Target mesh resolution in meters. This represents the characteristic
        edge length for mesh elements. The actual mesh area constraint is
        calculated as resolution squared.
    rift_name : str, optional
        Path to the file containing rift constraint geometry. If provided,
        these constraints will be incorporated into the mesh generation.
        Default is None.

    Returns
    -------
    md : object
        The input ISSM model object with updated mesh properties including:
        - mesh.x, mesh.y: Node coordinates
        - mesh.elements: Element connectivity matrix
        - mesh.segments: Boundary segment definitions
        - mesh.segmentmarkers: Boundary segment markers
        - mesh.numberofvertices: Total number of mesh vertices
        - mesh.numberofelements: Total number of mesh elements
        - mesh.vertexonboundary: Boolean array indicating boundary vertices
        - mesh.vertexconnectivity: Vertex-to-vertex connectivity
        - mesh.elementconnectivity: Element-to-element connectivity

    Raises
    ------
    IOError
        If the domain outline file or rift file (when specified) does not exist.
    RuntimeError
        If the Triangle Python wrappers are not installed and mesh creation fails.

    Warnings
    --------
    - Issues a warning if the existing mesh is not empty and will be overwritten.
    - Issues a warning if orphaned nodes are found and removed from the mesh.

    Notes
    -----
    This function requires the Triangle Python wrappers to be installed for
    mesh generation. If wrappers are not available, the function raises a RuntimeError.
    The function automatically handles orphaned nodes (nodes not belonging to
    any element) by removing them and updating the connectivity accordingly.

    Examples
    --------
    >>> import pyissm
    >>> md = pyissm.Model()
    >>> md = pyissm.tools.mesh.triangle(md, 'domain.exp', 1000.0)
    >>> md = pyissm.tools.mesh.triangle(md, 'domain.exp', 500.0, 'rifts.exp')
    """

    # Error checks
    ## Check if md mesh is empty
    if md.mesh.numberofelements:
        raise RuntimeError('md.mesh is not empty. Use md.mesh = pyissm.param.mesh.mesh2d() to reset the mesh.')

    ## Check if file(s) exist
    if not os.path.exists(domain_name):
        raise IOError(f"Domain outline file {domain_name} does not exist.")
    
    if rift_name is not None and not os.path.exists(rift_name):
        raise IOError(f"Rift file {rift_name} does not exist.")
    
    ## Check if wrappers are installed
    if not utils.wrappers.check_wrappers_installed():
        raise RuntimeError('pyissm.tools.mesh.triangle: Python wrappers not installed. Cannot create mesh.')


    # Calculate characteristic area. Resolution is node-oriented. A 1000 m resolution = 1000 * 1000 area
    area = resolution ** 2

    # Make the mesh using Triangle_python
    ## NOTE: Check for wrappers already done above
    elements, x, y, segments, segmentmarkers = utils.wrappers.Triangle(domain_name, rift_name, area)

    # Check that all created nodes belong to at least one element
    ## NOTE: Orphan node removal taken from $ISSM_DIR/src/m/mesh/triangle.m
    uniqueelements = np.sort(np.unique(elements))
    orphans = np.nonzero((~np.isin(range(1, len(x)), uniqueelements)).astype(int))[0]
    if len(orphans):
        warnings.warn(f'pyissm.tools.mesh.triangle: {len(orphans)} orphaned nodes found. These nodes do not belong to any element.\n'
                      'Removing orphaned nodes from the mesh.')
        
        for i in range(0, len(orphans)):
            print('WARNING: removing orphans')
            # Get rid of the orphan node i
            # Update x and y
            x = np.concatenate((x[0:(orphans[i] - i)], x[(orphans[i] - i + 1):]))
            y = np.concatenate((y[0:(orphans[i] - i)], y[(orphans[i] - i + 1):]))
            # Update elements
            pos = np.nonzero((elements > (orphans[i] - i)).flatten(order = 'F'))[0]
            elementstmp = elements.flatten(order = 'F')
            elementstmp[pos] -= 1
            elements = elementstmp.reshape(np.shape(elements), order = 'F')
            # Update segments
            pos1 = np.nonzero(segments[:,0] > (orphans[i] - i))[0]
            pos2 = np.nonzero(segments[:,1] > (orphans[i] - i))[0]
            segments[pos1, 0] -= 1
            segments[pos2, 1] -= 1
    
    # Assign to md structure
    md.mesh.x = x
    md.mesh.y = y
    md.mesh.elements = elements.astype(int)
    md.mesh.segments = segments.astype(int)
    md.mesh.segmentmarkers = segmentmarkers.astype(int)
    md.mesh.numberofvertices = len(x)
    md.mesh.numberofelements = len(elements)
    md.mesh.vertexonboundary = np.zeros(md.mesh.numberofvertices, int)
    md.mesh.vertexonboundary[md.mesh.segments[:,0:2] - 1] = 1

    ## Build connectivity arrays
    ### NOTE: Check for wrappers already done above
    md.mesh.vertexconnectivity = utils.wrappers.NodeConnectivity(md.mesh.elements, md.mesh.numberofvertices)
    md.mesh.elementconnectivity = utils.wrappers.ElementConnectivity(md.mesh.elements, md.mesh.vertexconnectivity)

    return md

def square_mesh(md,
                Lx,
                Ly,
                nx,
                ny):
    """
    Create a structured triangular mesh on a rectangular domain.
    
    This function generates a structured triangular mesh on a rectangular domain
    and populates the mesh fields of an ISSM model object. The mesh consists of
    triangular elements arranged in a regular grid pattern.

    Parameters
    ----------
    md : object
        ISSM model object whose mesh fields will be populated.
    Lx : float
        Length of the domain in the x-direction.
    Ly : float
        Length of the domain in the y-direction.
    nx : int
        Number of nodes in the x-direction.
    ny : int
        Number of nodes in the y-direction.
    
    Raises
    ------
    RuntimeError
        If Python wrappers are not installed.
    
    Warnings
    --------
    UserWarning
        If the model mesh is not empty and will be overwritten.
    
    Notes
    -----
    The function creates a structured triangular mesh by:
    1. Generating node coordinates on a regular grid
    2. Creating triangular elements by splitting each grid cell into two triangles
    3. Defining boundary segments around the domain perimeter
    4. Computing connectivity arrays for mesh topology
    The resulting mesh will have (nx-1)*(ny-1)*2 triangular elements and nx*ny nodes.

    Examples
    --------
    >>> import pyissm
    >>> md = pyissm.Model()
    >>> md = pyissm.tools.mesh.square_mesh(md, Lx=100.0, Ly=50.0, nx=11, ny=6)
    """  

    # Error checks
    ## Check if md mesh is empty
    if md.mesh.numberofelements:
        warnings.warn('md.mesh is not empty. Overwriting existing mesh.')

    ## Check if wrappers are installed
    if not utils.wrappers.check_wrappers_installed():
        raise RuntimeError('pyissm.tools.mesh.square_mesh: Python wrappers not installed. Cannot create mesh.')

    # Get number of elements and nodes
    n_elements = (nx - 1) * (ny - 1) * 2
    n_nodes = nx * ny

    # Initialize arrays
    index = np.zeros((n_elements, 3), int)
    x = np.zeros((n_nodes))
    y = np.zeros((n_nodes))

    # Create coordinates
    for n in range(0, nx):
        for m in range(0, ny):
            x[n * ny + m] = float(n)
            y[n * ny + m] = float(m)

    # Create index
    for n in range(0, nx - 1):
        for m in range(0, ny - 1):
            A = n * ny + (m + 1)
            B = A + 1
            C = (n + 1) * ny + (m + 1)
            D = C + 1
            index[n * (ny - 1) * 2 + 2 * m, :] = [A, C, B]
            index[n * (ny - 1) * 2 + 2 * (m + 1) - 1, :] = [B, C, D]

    # Scale x and y
    x = x / np.max(x) * Lx
    y = y / np.max(y) * Ly

    # Create segments
    segments = np.zeros((2 * (nx - 1) + 2 * (ny - 1), 3), int)
    segments[0:ny - 1, :] = np.vstack((np.arange(2, ny + 1), np.arange(1, ny), (2 * np.arange(1, ny) - 1))).T
    segments[ny - 1:2 * (ny - 1), :] = np.vstack((np.arange(ny * (nx - 1) + 1, nx * ny), np.arange(ny * (nx - 1) + 2, nx * ny + 1), 2 * np.arange((ny - 1) * (nx - 2) + 1, (nx - 1) * (ny - 1) + 1))).T
    segments[2 * (ny - 1):2 * (ny - 1) + (nx - 1), :] = np.vstack((np.arange(2 * ny, ny * nx + 1, ny), np.arange(ny, ny * (nx - 1) + 1, ny), np.arange(2 * (ny - 1), 2 * (nx - 1) * (ny - 1) + 1, 2 * (ny - 1)))).T
    segments[2 * (ny - 1) + (nx - 1):2 * (nx - 1) + 2 * (ny - 1), :] = np.vstack((np.arange(1, (nx - 2) * ny + 2, ny), np.arange(ny + 1, ny * (nx - 1) + 2, ny), np.arange(1, 2 * (nx - 2) * (ny - 1) + 2, 2 * (ny - 1)))).T

    # Populate md structure
    md.mesh.x = x
    md.mesh.y = y
    md.mesh.elements = index.astype(int)
    md.mesh.segments = segments.astype(int)
    md.mesh.numberofvertices = len(x)
    md.mesh.numberofelements = len(index)
    md.mesh.vertexonboundary = np.zeros(md.mesh.numberofvertices, int)
    md.mesh.vertexonboundary[md.mesh.segments[:,0:2] - 1] = 1

    ## Build connectivity arrays
    ### NOTE: Check for wrappers already done above
    md.mesh.vertexconnectivity = utils.wrappers.NodeConnectivity(md.mesh.elements, md.mesh.numberofvertices)
    md.mesh.elementconnectivity = utils.wrappers.ElementConnectivity(md.mesh.elements, md.mesh.vertexconnectivity)

    return md


def round_mesh(md,
               radius,
               resolution,
               exp_output_name = None,
               keep_exp = False):
    
    """
    Create a structured triangular mesh on a circular domain.
    
    This function generates a triangular mesh on a circular domain by first
    creating a domain outline file and then using the Triangle mesh generator.
    The mesh is created with approximately uniform resolution around the circle
    perimeter and moves the closest node to the origin.

    Parameters
    ----------
    md : object
        ISSM model object whose mesh fields will be populated.
    radius : float
        Radius of the circular domain in meters.
    resolution : float
        Target mesh resolution in meters. This represents the characteristic
        edge length for mesh elements around the circle perimeter.
    exp_output_name : str, optional
        Path for the output domain outline file (.exp format). If None,
        defaults to 'round_mesh.exp'. Default is None.
    keep_exp : bool, optional
        Whether to keep the temporary domain outline file after mesh creation.
        If False, the file is automatically deleted. Default is False.

    Returns
    -------
    md : object
        The input ISSM model object with updated mesh properties including:
        - mesh.x, mesh.y: Node coordinates
        - mesh.elements: Element connectivity matrix
        - mesh.segments: Boundary segment definitions
        - mesh.segmentmarkers: Boundary segment markers
        - mesh.numberofvertices: Total number of mesh vertices
        - mesh.numberofelements: Total number of mesh elements
        - mesh.vertexonboundary: Boolean array indicating boundary vertices
        - mesh.vertexconnectivity: Vertex-to-vertex connectivity
        - mesh.elementconnectivity: Element-to-element connectivity

    Raises
    ------
    IOError
        If the specified exp_output_name file already exists.
    RuntimeError
        If the Triangle Python wrappers are not installed.

    Warnings
    --------
    UserWarning
        If the model mesh is not empty and will be overwritten.

    Notes
    -----
    The function creates a circular mesh by:
    1. Generating points uniformly distributed around the circle perimeter
    2. Writing these points to a domain outline file (.exp format)
    3. Using the Triangle mesh generator to create the triangular mesh
    4. Moving the closest node to the origin (0,0) for convenience
    5. Optionally removing the temporary domain outline file

    The number of points on the circle perimeter is calculated based on the
    target resolution to ensure approximately uniform spacing.

    Examples
    --------
    >>> import pyissm
    >>> md = pyissm.Model()
    >>> md = pyissm.tools.mesh.round_mesh(md, radius=5000.0, resolution=500.0)
    >>> md = pyissm.tools.mesh.round_mesh(md, radius=1000.0, resolution=100.0, 
    ...                                   exp_output_name='circle.exp', keep_exp=True)
    """

    
    # Internal helper function
    def _round_sig_fig(x, n):
        nonzeros = np.where(x != 0)
        digits = np.ceil(np.log10(np.abs(x[nonzeros])))
        x[nonzeros] = x[nonzeros] / 10.**digits
        x[nonzeros] = np.round(x[nonzeros], decimals=n)
        x[nonzeros] = x[nonzeros] * 10.**digits
        return x
    
    ## ----------------------------------------------------------

    # Error checks
    ## NOTE: Existing mesh & wrapper installation checks handled in triangle()

    ## Check if file(s) exist. Do not overwrite existing files.
    if exp_output_name is not None and os.path.exists(exp_output_name):
        raise IOError(f"exp_output_name file {exp_output_name} already exists.")
    
    # Create the domain outline file
    if exp_output_name is None:
        exp_output_name = 'round_mesh.exp'

    # Construct the mesh
    ## Get number of points on the circle
    pointsonedge = int(np.floor((2. * np.pi * radius) / resolution) + 1)  # +1 to close the outline

    ## Calculate the Cartesian coordinates of the points
    theta = np.linspace(0., 2. * np.pi, pointsonedge)
    x_list = _round_sig_fig(radius * np.cos(theta), 12)
    y_list = _round_sig_fig(radius * np.sin(theta), 12)

    ## Create the contour dictionary
    contour = collections.OrderedDict()
    contour['x'] = x_list
    contour['y'] = y_list
    contour['density'] = 1.

    ## Write the contour to an exp file
    exp.exp_write(contour, exp_output_name)

    ## Create the mesh using triangle
    md = triangle(md, exp_output_name, resolution)

    ## Move closest node to (0,0)
    min_pos = np.argmin(np.sqrt(md.mesh.x**2 + md.mesh.y**2))
    md.mesh.x[min_pos] = 0.
    md.mesh.y[min_pos] = 0.

    # Remove temporary exp file if not required
    if not keep_exp:
        os.remove(exp_output_name)

    return md
    