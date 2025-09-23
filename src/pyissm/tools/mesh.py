"""
Tools for creating ISSM model meshes.
"""

import numpy as np
import warnings
import os
import collections

from .. import utils
from .. import param
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
    md.mesh = param.mesh.mesh2d()
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
    md.mesh = param.mesh.mesh2d()
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

def bamg(md, **kwargs):
    """
    Create a triangular mesh using the BAMG (Bidimensional Anisotropic Mesh Generator) algorithm.

    This function generates high-quality anisotropic triangular meshes for complex 
    geometries using the BAMG mesh generator. It supports various mesh constraints 
    including domain boundaries, holes, subdomains, rifts, and anisotropic metrics.

    Parameters
    ----------
    md : object
        ISSM model object whose mesh fields will be populated with the generated mesh.

    kwargs : dict, optional
        Additional keyword arguments to customize mesh generation. Supported options include:
        - anisomax : float, maximum anisotropy ratio allowed in the mesh (default= 1e30)
        - coeff : float, global mesh size coefficient multiplier (default= 1.0)
        - Crack : int, enable crack processing (0=disabled, 1=enabled) (default= 0)
        - cutoff : float, cutoff value for metric interpolation (default= 1e-5)
        - domain : str or list, path to domain outline file or list of domain contours (default= None)
        - err : float, target interpolation error for mesh adaptation (default= 0.01)
        - errg : float, target geometric error for mesh adaptation (default= 0.1)
        - field : array_like, field values for metric computation (default= empty array)
        - gradation : float, mesh size gradation parameter (default= 1.5)
        - Hessiantype : int, type of Hessian computation (0=P1, 1=P2) (default= 0)
        - hmax : float, maximum allowed edge length (default= 1e100)
        - hmin : float, minimum allowed edge length (default= 1e-100)
        - hmaxVertices : array_like, maximum edge lengths at specific vertices (default= empty array)
        - hminVertices : array_like, minimum edge lengths at specific vertices (default= empty array)
        - holes : str or list, path to holes file or list of hole contours (default= None)
        - hVertices : array_like, target edge lengths at specific vertices (default= empty array)
        - KeepVertices : int, keep vertices from previous mesh (0=no, 1=yes) (default= 1)
        - Markers : array_like, edge markers for boundary identification (default= None)
        - maxnbv : float, maximum number of vertices allowed (default= 1.0e6)
        - maxsubdiv : float, maximum number of edge subdivisions (default= 10.0)
        - metric : array_like, anisotropic metric tensor field (default= empty array)
        - Metrictype : int, type of metric (0=isotropic, 1=anisotropic) (default= 0)
        - nbjacobi : int, number of Jacobi smoothing iterations (default= 1)
        - nbsmooth : int, number of mesh smoothing iterations (default= 3)
        - NoBoundaryRefinement : int, disable boundary refinement for domain edges (0=allow, 1=disable) (default= 0)
        - NoBoundaryRefinementAllBoundaries : int, disable boundary refinement for all edges (0=allow, 1=disable) (default= 0)
        - omega : float, relaxation parameter for smoothing (default= 1.8)
        - power : float, power for metric computation (default= 1.0)
        - RequiredVertices : array_like, coordinates of vertices that must be included in the mesh (default= None)
        - rifts : str, path to rifts file for fracture modeling (default= None)
        - splitcorners : int, split corners in mesh generation (0=no, 1=yes) (default= 1)
        - subdomains : str or list, path to subdomains file or list of subdomain contours (default= None)
        - tol : float, tolerance for geometric operations (default= None)
        - toltip : float, tolerance for rift tip processing (default= None)
        - tracks : str or array_like, path to tracks file or track coordinates (default= None)
        - verbose : int, verbosity level (0=quiet, 1=verbose) (default= 1)
        - vertical : int, create 2D vertical mesh (0=standard 2D, 1=vertical) (default= 0)
        - 3dsurface : int, create 3D surface mesh (0=standard 2D, 1=3D surface) (default= 0)

    Returns
    -------
    md : object
        The input ISSM model object with updated mesh properties including:
        - mesh.x, mesh.y: Node coordinates
        - mesh.elements: Element connectivity matrix
        - mesh.edges: Edge connectivity matrix
        - mesh.segments: Boundary segment definitions
        - mesh.segmentmarkers: Boundary segment markers
        - mesh.numberofvertices: Total number of mesh vertices
        - mesh.numberofelements: Total number of mesh elements
        - mesh.numberofedges: Total number of mesh edges
        - mesh.vertexonboundary: Boolean array indicating boundary vertices
        - mesh.elementconnectivity: Element-to-element connectivity
        - private.bamg: BAMG-specific mesh and geometry data

    Raises
    ------
    IOError
        If specified input files (domain, holes, subdomains, rifts) do not exist.
    RuntimeError
        If mesh generation fails or if incompatible options are specified.
    TypeError
        If input arguments are of incorrect type.

    Notes
    -----
    This function is a comprehensive interface to the BAMG mesh generator, supporting
    complex geometries with multiple constraints. The mesh can be adapted based on
    metric fields for anisotropic meshing. Special handling is provided for rifts
    and fractures in ice sheet modeling applications.

    Examples
    --------
    >>> import pyissm
    >>> md = pyissm.Model()
    >>> # Basic domain meshing
    >>> md = pyissm.tools.mesh.bamg(md, domain='outline.exp', hmax=1000.0)
    >>> # Anisotropic meshing with metric
    >>> md = pyissm.tools.mesh.bamg(md, domain='outline.exp', 
    ...                             metric=metric_field, err=0.005)
    >>> # Mesh with holes and subdomains
    >>> md = pyissm.tools.mesh.bamg(md, domain='outline.exp', 
    ...                             holes='holes.exp', subdomains='regions.exp')
    """

    # Define default options
    defaults = {
        "anisomax": 1e30,
        "coeff": 1.0,
        "Crack": 0,
        "cutoff": 1e-5,
        "domain": None,
        "err": 0.01,
        "errg": 0.1,
        "field": np.empty((0, 1)),
        "gradation": 1.5,
        "Hessiantype": 0,
        "hmax": 1e100,
        "hmin": 1e-100,
        "hmaxVertices": np.empty((0, 1)),
        "hminVertices": np.empty((0, 1)),
        "holes": None,
        "hVertices": np.empty((0, 1)),
        "KeepVertices": 1,
        "Markers": None,
        "maxnbv": 1.0e6,
        "maxsubdiv": 10.0,
        "metric": np.empty((0, 1)),
        "Metrictype": 0,
        "nbjacobi": 1,
        "nbsmooth": 3,
        "NoBoundaryRefinement": 0,
        "NoBoundaryRefinementAllBoundaries": 0,
        "omega": 1.8,
        "power": 1.0,
        "RequiredVertices": None,
        "rifts": None,
        "splitcorners": 1,
        "subdomains": None,
        "tol": None,
        "toltip": None,
        "tracks": None,
        "verbose": 1,
        "vertical": 0,
        "3dsurface": 0
    }

    # Define helper functions
    def _bamg_geom(**kwargs):
        """
        Create a default BAMG geometry dictionary and update with user options.
        """
        geom = {
            'Vertices': np.empty((0, 3)),
            'Edges': np.empty((0, 3)),
            'TangentAtEdges': np.empty((0, 4)),
            'Corners': np.empty((0, 1)),
            'RequiredVertices': np.empty((0, 1)),
            'RequiredEdges': np.empty((0, 1)),
            'CrackedEdges': np.empty((0, 0)),
            'SubDomains': np.empty((0, 4)),
        }

        # Update defaults with user-specified options
        geom.update(kwargs)

        return geom

    def _bamg_mesh(**kwargs):
        """
        Create a default BAMG mesh dictionary and update with user options.
        """
        mesh = {
            'Vertices': np.empty((0, 3)),
            'Edges': np.empty((0, 3)),
            'Triangles': np.empty((0, 0)),
            'IssmEdges': np.empty((0, 0)),
            'IssmSegments': np.empty((0, 0)),
            'VerticesOnGeomVertex': np.empty((0, 0)),
            'VerticesOnGeomEdge': np.empty((0, 0)),
            'EdgesOnGeomEdge': np.empty((0, 0)),
            'SubDomains': np.empty((0, 4)),
            'SubDomainsFromGeom': np.empty((0, 0)),
            'ElementConnectivity': np.empty((0, 0)),
            'NodalConnectivity': np.empty((0, 0)),
            'NodalElementConnectivity': np.empty((0, 0)),
            'CrackedVertices': np.empty((0, 0)),
            'CrackedEdges': np.empty((0, 0))
        }

        # Update defaults with user-specified options
        mesh.update(kwargs)

        return mesh

    def _load_spatial_components(component):
        """
        Load spatial components from file or directly from input.
        """

        if component is not None:

            ## Check the file exists if a filename is provided
            if isinstance(component, str):
                if not os.path.exists(component):
                    raise IOError(f"BAMG spatial component file {component} does not exist.")
                return exp.exp_read(component)
            
            ## If a list or dict is provided, it must be a list of dictionaries
            elif isinstance(component, list):
                if len(component):
                    if all(isinstance(c, (dict, collections.OrderedDict)) for c in component):
                        return component
                    else:
                        raise Exception("pyissm.tools.mesh.bamg: if a spatial component is a list, its elements must be of type dict or OrderedDict")
                else:
                    return component
            
            ## Single contour provided as a dict, return list with one element
            elif isinstance(component, (dict, collections.OrderedDict)):
                return [component]
            else:
                raise Exception("pyissm.tools.mesh.bamg: spatial components must be a filename (str), a list of contours (list of dicts), or a single contour (dict).")
        return []
    
    def _process_spatial_component(component, domain, ref_counter, count, is_hole = False, is_subdomain = False):
        """
        Process a spatial component (domain, holes, subdomains) for BAMG.
        """

        ## For each contour in the component
        for i in range(len(component)):

            ### Check the contour is closed
            if component[i]['x'][0] != component[i]['x'][-1] or component[i]['y'][0] != component[i]['y'][-1]:
                raise Exception("pyissm.tools.mesh.bamg: each contour must be closed (first and last points must be identical).")

            ### Check all holes and subdomain contours are INSIDE the principle domain
            ## NOTE: This check differs from existing bamg.m and bamg.py which only checks the domain contours (not holes)
            if is_hole or is_subdomain:
                for c in component:
                    flags = utils.wrappers.ContourToNodes(c['x'], c['y'], [domain[0]], 0)[0]
                    if np.any(np.logical_not(flags)):
                        raise Exception("pyissm.tools.mesh.bamg: all contours in 'holes' and 'subdomains' must be inside the first contour of 'domain'.")

            ### Check the orientation of the contour
            nods = component[i]['nods'] - 1
            test = np.sum((component[i]['x'][1:nods + 1] - component[i]['x'][0:nods]) * (component[i]['y'][1:nods + 1] + component[i]['y'][0:nods]))

            if (is_hole and test < 0) or (is_subdomain and test > 0):
                ## TODO: Why is subdomain orientation opposite to holes?
                print("At least one contour in 'holes' or 'subdomains' was not correctly oriented and has been re-oriented")
                component[i]['x'] = np.flipud(component[i]['x'])
                component[i]['y'] = np.flipud(component[i]['y'])

            ## If processing domain, add info to bamg_geom
            if not is_hole and not is_subdomain:
                ### Flag how many edges we have so far
                edge_length = len(bamg_geom['Edges'])

                ### Add all points to bamg_geom
                bamg_geom['Vertices'] = np.vstack((bamg_geom['Vertices'], np.vstack((component[i]['x'][0:nods], component[i]['y'][0:nods], np.ones((nods)))).T))
                bamg_geom['Edges'] = np.vstack((bamg_geom['Edges'], np.vstack((np.arange(count + 1, count + nods + 1), np.hstack((np.arange(count + 2, count + nods + 1), count + 1)), 1. * np.ones((nods)))).T))

                new_edge_length = len(bamg_geom['Edges'])
                edges_required = np.asarray(range((edge_length + 1), (new_edge_length + 1)))  # NOTE: Upper bound of range is non-inclusive (compare to src/m/mesh/bamg.m)
                if i > 0:
                    bamg_geom['SubDomains'] = np.vstack((bamg_geom['SubDomains'], [2, count + 1, 1, -ref_counter]))
                    ref_counter = ref_counter + 1
                else:
                    bamg_geom['SubDomains'] = np.vstack((bamg_geom['SubDomains'], [2, count + 1, 1, 0]))

            ## If processing holes, add info to bamg_geom
            ## TODO: Why negative ref_counter for holes?
            if is_hole:
                bamg_geom['Vertices'] = np.vstack((bamg_geom['Vertices'], np.vstack((holes[i]['x'][0:nods], holes[i]['y'][0:nods], np.ones((nods)))).T))
                bamg_geom['Edges'] = np.vstack((bamg_geom['Edges'], np.vstack((np.arange(count + 1, count + nods + 1), np.hstack((np.arange(count + 2, count + nods + 1), count + 1)), 1. * np.ones((nods)))).T))
                bamg_geom['SubDomains'] = np.vstack((bamg_geom['SubDomains'], [2, count + 1, 1, -ref_counter]))

                ref_counter = ref_counter + 1
            
            ## If processing subdomains, add info to bamg_geom
            ## TODO: Why positive ref_counter for subdomains?
            if is_subdomain:
                bamg_geom['Vertices'] = np.vstack((bamg_geom['Vertices'], np.vstack((subdomains[i]['x'][0:nods], subdomains[i]['y'][0:nods], np.ones((nods)))).T))
                bamg_geom['Edges'] = np.vstack((bamg_geom['Edges'], np.vstack((np.arange(count + 1, count + nods + 1), np.hstack((np.arange(count + 2, count + nods + 1), count + 1)), 1. * np.ones((nods)))).T))
                bamg_geom['SubDomains'] = np.vstack((bamg_geom['SubDomains'], [2, count + 1, 1, ref_counter]))

                ref_counter = ref_counter + 1

            # Update counter
            count += nods

        if not (is_hole or is_subdomain):
            return component, ref_counter, count, bamg_geom, edges_required
        else:
            return component, ref_counter, count, bamg_geom
    
    def _seg_intersect(seg1, seg2):
        """
        Check if 2 segments intersect.
        NOTE: Taken directly from $ISSM_DIR/src/m/geometry/SegIntersect.py
        """

        bval = 1

        xA = seg1[0, 0]
        yA = seg1[0, 1]
        xB = seg1[1, 0]
        yB = seg1[1, 1]
        xC = seg2[0, 0]
        yC = seg2[0, 1]
        xD = seg2[1, 0]
        yD = seg2[1, 1]

        O2A = np.array([xA, yA]) - np.array([xD / 2. + xC / 2., yD / 2. + yC / 2.])
        O2B = np.array([xB, yB]) - np.array([xD / 2. + xC / 2., yD / 2. + yC / 2.])
        O1C = np.array([xC, yC]) - np.array([xA / 2. + xB / 2., yB / 2. + yA / 2.])
        O1D = np.array([xD, yD]) - np.array([xA / 2. + xB / 2., yB / 2. + yA / 2.])

        n1 = np.array([yA - yB, xB - xA])  #normal vector to segA
        n2 = np.array([yC - yD, xD - xC])  #normal vector to segB

        test1 = np.dot(n2, O2A)
        test2 = np.dot(n2, O2B)

        if test1 * test2 > 0:
            bval = 0
            return bval

        test3 = np.dot(n1, O1C)
        test4 = np.dot(n1, O1D)

        if test3 * test4 > 0:
            bval = 0
            return bval

        #if colinear
        if test1 * test2 == 0 and test3 * test4 == 0 and np.linalg.det(np.hstack((n1.reshape((-1, )), n2.reshape(-1, )))) == 0:

            #projection on the axis O1O2
            O2O1 = np.array([xA / 2. + xB / 2., yB / 2. + yA / 2.]) - np.array([xD / 2. + xC / 2., yD / 2. + yC / 2.])
            O1A = np.dot(O2O1, (O2A - O2O1))
            O1B = np.dot(O2O1, (O2B - O2O1))
            O1C = np.dot(O2O1, O1C)
            O1D = np.dot(O2O1, O1D)

        #test if one point is included in the other segment (-> bval = 1)
            if (O1C - O1A) * (O1D - O1A) < 0:
                bval = 1
                return bval
            if (O1C - O1B) * (O1D - O1B) < 0:
                bval = 1
                return bval
            if (O1A - O1C) * (O1B - O1C) < 0:
                bval = 1
                return bval
            if (O1A - O1D) * (O1B - O1D) < 0:
                bval = 1
                return bval

        #test if the 2 segments have the same middle (-> bval = 1)
            if O2O1 == 0:
                bval = 1
                return bval

        #else
            bval = 0
            return bval

        return bval

    def _process_rifts(rift_file):
        """
        Process rifts for BAMG.
        """
        
        # Error checks
        if not isinstance(rift_file, (str)):
            raise TypeError("pyissm.tools.mesh.bamg: rifts must be a filename (str).")
        if not os.path.exists(rift_file):
            raise IOError(f"pyissm.tools.mesh.bamg: rift file {rift_file} does not exist.")

        # Read rift file
        rift = exp.exp_read(rift_file)

        # Process each rift
        for i in range(len(rift)):

            ## Check whether all points of the rift are inside the domain
            flags = utils.wrappers.ContourToNodes(rift[i]['x'], rift[i]['y'], [domain[0]], 0)[0]
            if np.all(np.logical_not(flags)):
                raise RuntimeError("pyissm.tools.mesh.bamg: one rift has all its points outside of the domain outline")
            
            ## Check if rift tip is outside of the domain
            elif np.any(np.logical_not(flags)):
                # We have LOTS of work to do
                print('Rift tip outside of or on the domain has been detected and is being processed...')

                ## Check that only one point is outside (for now)
                if np.sum(np.logical_not(flags).astype(int)) != 1:
                    raise RuntimeError("pyissm.tools.mesh.bamg: only one point outside of the domain is supported at this time")

                ## Move tip outside to the first position
                if not flags[0]:
                    # OK, first point is outside (do nothing),
                    pass
                elif not flags[-1]:
                    rift[i]['x'] = np.flipud(rift[i]['x'])
                    rift[i]['y'] = np.flipud(rift[i]['y'])
                else:
                    raise RuntimeError('pyissm.tools.mesh.bamg: only a rift tip can be outside of the domain')

                # Get coordinate of intersection point
                x1 = rift[i]['x'][0]
                y1 = rift[i]['y'][0]
                x2 = rift[i]['x'][1]
                y2 = rift[i]['y'][1]
                for j in range(0, np.size(domain[0]['x']) - 1):
                    if _seg_intersect(np.array([[x1, y1], [x2, y2]]), np.array([[domain[0]['x'][j], domain[0]['y'][j]], [domain[0]['x'][j + 1], domain[0]['y'][j + 1]]])):

                        # Get position of the two nodes of the edge in domain
                        i1 = j
                        i2 = j + 1

                        # Rift is crossing edge [i1, i2] of the domain
                        # Get coordinate of intersection point (http://mathworld.wolfram.com/Line-LineIntersection.html)
                        x3 = domain[0]['x'][i1]
                        y3 = domain[0]['y'][i1]
                        x4 = domain[0]['x'][i2]
                        y4 = domain[0]['y'][i2]
                        x = np.linalg.det(np.array([[np.linalg.det(np.array([[x1, y1], [x2, y2]])), x1 - x2], [np.linalg.det(np.array([[x3, y3], [x4, y4]])), x3 - x4]])) / np.linalg.det(np.array([[x1 - x2, y1 - y2], [x3 - x4, y3 - y4]]))
                        y = np.linalg.det(np.array([[np.linalg.det(np.array([[x1, y1], [x2, y2]])), y1 - y2], [np.linalg.det(np.array([[x3, y3], [x4, y4]])), y3 - y4]])) / np.linalg.det(np.array([[x1 - x2, y1 - y2], [x3 - x4, y3 - y4]]))

                        segdis = np.sqrt((x4 - x3)**2 + (y4 - y3)**2)
                        tipdis = np.array([np.sqrt((x - x3)**2 + (y - y3)**2), np.sqrt((x - x4)**2 + (y - y4)**2)])

                        if np.min(tipdis) / segdis < options['toltip']:
                            print('moving tip-domain intersection point')

                            # Get position of the closer point
                            if tipdis[0] > tipdis[1]:
                                pos = i2
                            else:
                                pos = i1

                            # This point is only in Vertices (number pos).
                            # OK, now we can add our own rift
                            nods = rift[i]['nods'] - 1
                            bamg_geom['Vertices'] = np.vstack((bamg_geom['Vertices'], np.hstack((rift[i]['x'][1:].reshape(-1, ), rift[i]['y'][1:].reshape(-1, ), np.ones((nods, 1))))))
                            bamg_geom['Edges'] = np.vstack((
                                bamg_geom['Edges'],
                                np.array([[pos, count + 1, (1 + i)]]),
                                np.hstack((np.arange(count + 1, count + nods).reshape(-1, ), np.arange(count + 2, count + nods + 1).reshape(-1, ), (1 + i) * np.ones((nods - 1, 1))))
                            ))
                            count += nods
                            break
                        else:
                            # Add intersection point to Vertices
                            bamg_geom['Vertices'] = np.vstack((bamg_geom['Vertices'],
                                np.array([[x, y, 1]])
                            ))
                            count += 1

                            # Decompose the crossing edge into 2 subedges
                            pos = np.nonzero(np.logical_and(bamg_geom['Edges'][:, 0] == i1, bamg_geom['Edges'][:, 1] == i2))[0]
                            if not pos:
                                raise RuntimeError('pyissm.tools.mesh.bamg: a problem occurred...')
                            bamg_geom['Edges'] = np.vstack((
                                bamg_geom['Edges'][0:pos - 1, :],
                                np.array([[
                                    bamg_geom['Edges'][pos, 0],
                                    count,
                                    bamg_geom['Edges'][pos, 2]
                                ]]),
                                np.array([[
                                    count,
                                    bamg_geom['Edges'][pos, 1],
                                    bamg_geom['Edges'][pos, 2]
                                ]]),
                                bamg_geom['Edges'][pos + 1:, :]
                            ))

                            # OK, now we can add our own rift
                            nods = rift[i]['nods'] - 1
                            bamg_geom['Vertices'] = np.vstack((bamg_geom['Vertices'],
                                np.hstack((
                                    rift[i]['x'][1:].reshape(-1, ),
                                    rift[i]['y'][1:].reshape(-1, ),
                                    np.ones((nods, 1))
                                ))
                            ))
                            bamg_geom['Edges'] = np.vstack((
                                bamg_geom['Edges'],
                                np.array([[count, count + 1, 2]]),
                                np.hstack((
                                    np.arange(count + 1, count + nods).reshape(-1, ),
                                    np.arange(count + 2, count + nods + 1).reshape(-1, ),
                                    (1 + i) * np.ones((nods - 1, 1))
                                ))
                            ))
                            count += nods
                            break
            else:
                nods = rift[i]['nods'] - 1
                bamg_geom['Vertices'] = np.vstack((
                    bamg_geom['Vertices'],
                    np.hstack((
                        rift[i]['x'][:],
                        rift[i]['y'][:],
                        np.ones((nods + 1, 1))
                    ))
                ))
                bamg_geom['Edges'] = np.vstack((
                    bamg_geom['Edges'],
                    np.hstack((
                        np.arange(count + 1, count + nods).reshape(-1, ),
                        np.arange(count + 2, count + nods + 1).reshape(-1, ),
                        i * np.ones((nods, 1))
                    ))
                ))
                count += (nods + 1)

        return count, bamg_geom
    
    def _process_tracks(track):
        """
        Process tracks for BAMG.
        """
                
        # Read tracks
        if all(isinstance(track, str)):
            track = exp.expread(track)
            track = np.hstack((track.x.reshape(-1, ), track.y.reshape(-1, )))
        else:
            track = float(track)

        if np.size(track, axis=1) == 2:
            track = np.hstack((track, 3. * np.ones((np.size(track, axis=0), 1))))

        # Only keep those inside
        flags = utils.wrappers.ContourToNodes(track[:, 0], track[:, 1], [domain[0]], 0)[0]
        track = track[np.nonzero(flags), :]

        # Add all points to bamg_geometry
        nods = np.size(track, axis=0)
        bamg_geom['Vertices'] = np.vstack((bamg_geom['Vertices'], track))
        bamg_geom['Edges'] = np.vstack((
            bamg_geom['Edges'],
            np.hstack((
                np.arange(count + 1, count + nods).reshape(-1, ),
                np.arange(count + 2, count + nods + 1).reshape(-1, ),
                3. * np.ones((nods - 1, 1))
            ))
        ))

        # Update counter
        count += nods

        return count, bamg_geom
    
    def _process_required_vertices(required_vertices):
        """
        Process required vertices for BAMG.
        """

        if np.size(required_vertices, axis=1) == 2:
            required_vertices = np.hstack((required_vertices, 4. * np.ones((np.size(required_vertices, axis=0), 1))))

        # Only keep those inside
        flags = utils.wrappers.ContourToNodes(required_vertices[:, 0], required_vertices[:, 1], [domain[0]], 0)[0]
        required_vertices = required_vertices[np.nonzero(flags)[0], :]

        # Add all points to bamg_geom
        nods = np.size(required_vertices, axis=0)
        bamg_geom['Vertices'] = np.vstack((bamg_geom['Vertices'], required_vertices))

        # Update counter
        count += nods

        return count, bamg_geom
    
    # ---------------------------------------------------------------
    
    # Update defaults with user-specified options
    options = collections.OrderedDict(defaults)
    options.update(kwargs)

    # Initialize geometry and mesh structures
    bamg_geom = _bamg_geom()
    bamg_mesh = _bamg_mesh()

    # Initialise counters
    subdomain_ref = 1
    hole_ref = 1

    # Build BAMG mesh from domain, holes, subdomains
    if options['domain'] is not None:
        domain_file = options['domain']
        domain = _load_spatial_components(domain_file)

        holes = []
        if options['holes'] is not None:
            hole_file = options['holes']
            holes = _load_spatial_components(hole_file)

        subdomains = []
        if options['subdomains'] is not None:
            subdomain_file = options['subdomains']
            subdomains = _load_spatial_components(subdomain_file)

        # Build geometry
        count = 0

        ## Process domain
        domain, subdomain_ref, count, bamg_geom, edges_required = _process_spatial_component(domain, domain, subdomain_ref, count, is_hole = False, is_subdomain = False)

        ## Process holes
        holes, hole_ref, count, bamg_geom = _process_spatial_component(holes, domain, hole_ref, count, is_hole = True, is_subdomain = False)

        ## Process subdomains
        subdomains, subdomain_ref, count, bamg_geom = _process_spatial_component(subdomains, domain, subdomain_ref, count, is_hole = False, is_subdomain = True)

        # Process vertical options
        if options['vertical'] == 1:
            if np.size(options['Markers']) != np.size(bamg_geom['Edges'], 0):
                edges_size = np.size(bamg_geom['Edges'], 0)
                raise RuntimeError(f'for 2d vertical mesh, \'Markers\' option is required, and should be of size {edges_size}')
        
        if np.size(options['Markers']) == np.size(bamg_geom['Edges'], 0):
            bamg_geom['Edges'][:, 2] = options['Markers']

        # Process rifts
        if options['rifts'] is not None:
            count, bamg_geom = _process_rifts(options['rifts'])

        # Process tracks
        if options['tracks'] is not None:
            count, bamg_geom = _process_tracks(options['tracks'])

        # Proces required vertices
        if options['RequiredVertices'] is not None:
            count, bamg_geom = _process_required_vertices(options['RequiredVertices'])

        # Process RequiredEdges
        if options['NoBoundaryRefinement'] == 1:
            bamg_geom['RequiredEdges'] = edges_required
        elif options['NoBoundaryRefinementAllBoundaries'] == 1:
                bamg_geom['RequiredEdges'] = np.arange(1, bamg_geom['Edges'].shape[0]).T
    
    # If a geometry is already provided, use it
    elif isinstance(md.private.bamg, dict) and 'geometry' in md.private.bamg:
        bamg_geom = _bamg_geom(**md.private.bamg['geometry'])
    else:
        # Do nothing...
        pass

    # If domain is not specified, check for existing mesh
    if md.mesh.numberofvertices and md.mesh.element_type() == 'Tria':
        ## If there is an existing BAMG mesh, use it
        if isinstance(md.private.bamg, dict) and 'mesh' in md.private.bamg:
            bamg_mesh = _bamg_mesh(**md.private.bamg['mesh'])
        else:
            ## If there is an existing non-BAMG mesh, convert it
            bamg_mesh['Vertices'] = np.vstack((
                md.mesh.x,
                md.mesh.y,
                np.ones((md.mesh.numberofvertices))
            )).T
            bamg_mesh['Triangles'] = np.hstack((md.mesh.elements, np.ones((md.mesh.numberofelements, 1))))

        ## If there are rifts in the model, raise an error (not supported yet)
        if isinstance(md.rifts.riftstruct, dict):
            raise TypeError('pyissm.tools.mesh.bamg: rifts not supported yet. Do meshprocessrift after bamg.')
    
    # Call the BAMG mesher
    bamg_mesh_out, bamg_geom_out = utils.wrappers.BamgMesher(bamg_mesh, bamg_geom, options)

    # Populate md structure
    if options['vertical'] == 1:
        ## Create 2D vertical mesh
        md.mesh = param.mesh.mesh2dvertical()
        md.mesh.x = bamg_mesh_out['Vertices'][:, 0].copy()
        md.mesh.y = bamg_mesh_out['Vertices'][:, 1].copy()
        md.mesh.elements = bamg_mesh_out['Triangles'][:, 0:3].astype(int)
        md.mesh.edges = bamg_mesh_out['IssmEdges'].astype(int)
        md.mesh.segments = bamg_mesh_out['IssmSegments'][:, 0:3].astype(int)
        md.mesh.segmentmarkers = bamg_mesh_out['IssmSegments'][:, 3].astype(int)

        md.mesh.numberofelements = np.size(md.mesh.elements, axis=0)
        md.mesh.numberofvertices = np.size(md.mesh.x)
        md.mesh.numberofedges = np.size(md.mesh.edges, axis=0)
        md.mesh.vertexonboundary = np.zeros(md.mesh.numberofvertices, int)
        md.mesh.vertexonboundary[md.mesh.segments[:, 0:2] - 1] = 1

    elif options['3dsurface'] == 1:
        md.mesh = param.mesh.mesh3dsurface()
        md.mesh.x = bamg_mesh_out['Vertices'][:, 0].copy()
        md.mesh.y = bamg_mesh_out['Vertices'][:, 1].copy()
        md.mesh.z = md.mesh.x
        md.mesh.z[:] = 0
        md.mesh.elements = bamg_mesh_out['Triangles'][:, 0:3].astype(int)
        md.mesh.edges = bamg_mesh_out['IssmEdges'].astype(int)
        md.mesh.segments = bamg_mesh_out['IssmSegments'][:, 0:3].astype(int)
        md.mesh.segmentmarkers = bamg_mesh_out['IssmSegments'][:, 3].astype(int)

        # Fill in rest of fields
        md.mesh.numberofelements = np.size(md.mesh.elements, axis=0)
        md.mesh.numberofvertices = np.size(md.mesh.x)
        md.mesh.numberofedges = np.size(md.mesh.edges, axis=0)
        md.mesh.vertexonboundary = np.zeros(md.mesh.numberofvertices, int)
        md.mesh.vertexonboundary[md.mesh.segments[:, 0:2] - 1] = 1

    else:
        md.mesh = param.mesh.mesh2d()
        md.mesh.x = bamg_mesh_out['Vertices'][:, 0].copy()
        md.mesh.y = bamg_mesh_out['Vertices'][:, 1].copy()
        md.mesh.elements = bamg_mesh_out['Triangles'][:, 0:3].astype(int)
        md.mesh.edges = bamg_mesh_out['IssmEdges'].astype(int)
        md.mesh.segments = bamg_mesh_out['IssmSegments'][:, 0:3].astype(int)
        md.mesh.segmentmarkers = bamg_mesh_out['IssmSegments'][:, 3].astype(int)

        # Fill in rest of fields
        md.mesh.numberofelements = np.size(md.mesh.elements, axis=0)
        md.mesh.numberofvertices = np.size(md.mesh.x)
        md.mesh.numberofedges = np.size(md.mesh.edges, axis=0)
        md.mesh.vertexonboundary = np.zeros(md.mesh.numberofvertices, int)
        md.mesh.vertexonboundary[md.mesh.segments[:, 0:2] - 1] = 1
    

    # BAMG private fields
    md.private.bamg = collections.OrderedDict()
    md.private.bamg['mesh'] = _bamg_mesh(**bamg_mesh_out)
    md.private.bamg['geometry'] = _bamg_geom(**bamg_geom_out)
    md.mesh.elementconnectivity = md.private.bamg['mesh']['ElementConnectivity']
    md.mesh.elementconnectivity[np.nonzero(np.isnan(md.mesh.elementconnectivity))] = 0
    md.mesh.elementconnectivity = md.mesh.elementconnectivity.astype(int)

    # Check for orphan vertices
    if np.any(np.logical_not(np.isin(np.arange(1, md.mesh.numberofvertices + 1), md.mesh.elements.flat))):
        raise RuntimeError('Output mesh has orphans. Check your Domain and/or RequiredVertices')

    return md

def bamgflowband(md,
                 x,
                 surf,
                 base,
                 **kwargs):
    
    """
    Create a flowband mesh using BAMG (Bidimensional Anisotropic Mesh Generator).

    This function generates a triangular mesh for a flowband (vertical 2D slice) using 
    the BAMG mesh generator. The flowband is defined by surface and base profiles 
    along a specified coordinate path, creating a vertical mesh suitable for ice flow 
    modeling in the vertical plane.

    Parameters
    ----------
    md : object
        ISSM model object whose mesh fields will be populated with the generated flowband mesh.
    x : array_like
        1D array of coordinates along the flowband path. These represent the horizontal
        positions where the surface and base elevations are defined.
    surf : array_like
        1D array of surface elevations corresponding to each x coordinate. Must have
        the same length as x.
    base : array_like
        1D array of base (bed) elevations corresponding to each x coordinate. Must have
        the same length as x.
    **kwargs : dict, optional
        Additional keyword arguments passed to the bamg function. See bamg() documentation
        for supported options.

    Returns
    -------
    md : object
        A new ISSM model object with a vertical 2D mesh populated, including:
        - mesh.x, mesh.y: Node coordinates in the flowband coordinate system
        - mesh.elements: Element connectivity matrix for triangular elements
        - mesh.edges: Edge connectivity matrix
        - mesh.segments: Boundary segment definitions with markers
        - mesh.segmentmarkers: Boundary segment markers (1=base, 2=right, 3=surface, 4=left)
        - mesh.numberofvertices: Total number of mesh vertices
        - mesh.numberofelements: Total number of mesh elements
        - mesh.numberofedges: Total number of mesh edges
        - mesh.vertexonboundary: Boolean array indicating boundary vertices
        - mesh.vertexonbase: Boolean array indicating vertices on the base boundary
        - mesh.vertexonsurface: Boolean array indicating vertices on the surface boundary
        - mesh.elementconnectivity: Element-to-element connectivity

    Raises
    ------
    ValueError
        If x, surf, and base arrays do not have the same length.
    RuntimeError
        If BAMG mesh generation fails or if incompatible meshing options are specified.

    Notes
    -----
    This function creates a vertical 2D mesh by:
    1. Constructing a closed domain from the surface and base profiles
    2. Assigning boundary markers: 1=base, 2=right side, 3=surface, 4=left side
    3. Calling the BAMG mesh generator with vertical=1 option (to convert to 2D vertical mesh)
    4. Post-processing to identify vertices on base and surface boundaries

    The resulting mesh is suitable for flowband modeling where ice flow is assumed
    to be primarily in the vertical plane defined by the x-coordinate path.

    Examples
    --------
    >>> import numpy as np
    >>> import pyissm
    >>> md = pyissm.Model()
    >>> x = np.arange(1, 3001, 100).T
    >>> h = np.linspace(1000, 300, np.size(x)).T
    >>> b = -917. / 1023. * h
    >>> md = pyissm.tools.mesh.bamgflowband(md, x = x, surf = b + h, base = b, hmax = 80.)
    """

    # Create domain structure
    domain = collections.OrderedDict()
    domain['x'] = np.concatenate((x, np.flipud(x), [x[0]]))
    domain['y'] = np.concatenate((base, np.flipud(surf), [base[0]]))
    domain['nods'] = np.size(domain['x'])

    # Create markers (base, right side, top surface, left side)
    m = np.ones((np.size(domain['x']) - 1, ))
    m[np.size(x) - 1] = 2
    m[np.size(x):2 * np.size(x) - 1] = 3
    m[2 * np.size(x) - 1] = 4

    # Call bamg
    md = core.Model()
    md = bamg(md, domain = [domain], Markers = m, vertical = 1, **kwargs)

    # Deal with vertices on bed
    ## NOTE: vertexonbase and vertexonsurface used to be set using vertexflags() defined in mesh2dvertical.py
    ## Here, we just do this inline because it's only used here and it's simpler this way.
    md.mesh.vertexonbase = np.zeros((md.mesh.numberofvertices, ))
    base_segments = md.mesh.segments[np.where(md.mesh.segmentmarkers == 1), 0:2] - 1
    md.mesh.vertexonbase[base_segments] = 1

    md.mesh.vertexonsurface = np.zeros((md.mesh.numberofvertices, ))
    surface_segments = md.mesh.segments[np.where(md.mesh.segmentmarkers == 3), 0:2] - 1
    md.mesh.vertexonsurface[surface_segments] = 1

    return md

    

