"""
Tools working with *.exp files for ISSM model domains and contours.
"""

import numpy as np
import collections
import os

def exp_write(contours, filename):
    """
    Write contours to an exp file.

    This function writes contour data to a file in *.exp format. The function can handle
    both single contours and lists of contours. Each contour should be a dictionary
    containing 'x' and 'y' coordinate data, and optionally 'name' and 'density' fields.

    Parameters
    ----------
    contours : dict or list of dict
        Contour data to write. If a dictionary, represents a single contour.
        If a list, represents multiple contours. Each contour dictionary should contain:
        - 'x' : array-like or scalar
            X coordinates of the contour points
        - 'y' : array-like or scalar  
            Y coordinates of the contour points
        - 'name' : str, optional
            Name of the contour. If not provided, filename is used
        - 'density' : int, optional
            Density value for the contour. Default is 1.0
    filename : str
        Path to the output exp file

    Raises
    ------
    RuntimeError
        If X and Y coordinates are not of identical size

    Notes
    -----
    The exp file format includes headers with contour names, point counts,
    density values, and coordinate data formatted to 10 decimal places.

    Examples
    --------
    >>> contour = {'x': [0, 1, 2], 'y': [0, 1, 0], 'name': 'triangle'}
    >>> exp_write(contour, 'output.exp')
    >>> contours = [
    ...     {'x': [0, 1], 'y': [0, 1], 'density': 2},
    ...     {'x': [2, 3], 'y': [2, 3], 'name': 'line2'}
    ... ]
    >>> exp_write(contours, 'multiple_contours.exp')
    """
    # Internal helper functions
    def _write_geom_list(contour, fid, filename):
        
        # Error checks
        if len(contour['x']) != len(contour['y']):
            raise RuntimeError('pyissm.tools.exp.exp_write: X and Y coordinates must be of identical size.')
        
        # Write contour to file
        ## If contour has a name, use it. Otherwise use the filename
        if 'name' in contour:
            fid.write('{}{}\n'.format('## Name:', contour['name']))
        else:
            fid.write('{}{}\n'.format('## Name:', filename))

        ## Write header information for the contour
        fid.write('{}\n'.format('## Icon:0'))
        fid.write('{}\n'.format('# Points Count Value'))
        
        ## Write point count and density value
        if 'density' in contour and isinstance(contour['density'], (int)):
            if isinstance(contour['density'], int):
                fid.write('{} {}\n'.format(np.size(contour['x']), contour['density']))
        else:
            ### Use default density of 1.0 if no density specified (or it's not an integer)
            fid.write('{} {}\n'.format(np.size(contour['x']), 1.))
        
        ## Write coordinate data header
        fid.write('{}\n'.format('# X pos Y pos'))
        
        ## Write each coordinate pair
        for x, y in zip(contour['x'], contour['y']):
            fid.write('%10.10f %10.10f\n' % (x, y))
        
        ## Add blank line after contour
        fid.write('\n')

    def _write_geom(contour, fid, filename):

        # Error checks
        if len(contour['x']) != len(contour['y']):
            raise RuntimeError('pyissm.tools.exp.exp_write: X and Y coordinates must be of identical size.')
        
        # Write contour to file
        ## If contour has a name, use it. Otherwise use the filename        
        if 'name' in contour:
            fid.write('{}{}\n'.format('## Name:', contour['name']))
        else:
            fid.write('{}{}\n'.format('## Name:', filename))

        ## Write header information for the contour
        fid.write('{}\n'.format('## Icon:0'))
        fid.write('{}\n'.format('# Points Count Value'))

        ## Write point count and density value
        if 'density' in contour and isinstance(contour['density'], (int)):
            if isinstance(contour['density'], int):
                fid.write('{} {}\n'.format(1, contour['density']))
        else:
            ### Use default density of 1.0 if no density specified (or it's not an integer)
            fid.write('{} {}\n'.format(1, 1.))

        ## Write coordinate data header
        fid.write('{}\n'.format('# X pos Y pos'))

        ## Write coordinate pairs
        fid.write('%10.10f %10.10f\n' % (contour['x'], contour['y']))

        ## Add blank line after contour
        fid.write('\n')

    ## ----------------------------------------------------------
    
    # Open the file for writing
    fid = open(filename, 'w')
    
    # If contours is a list, loop over several contours
    if isinstance(contours, list):
        for contour in contours:
            ## If contour is an array, loop on indexes
            if isinstance(contour['x'], (list, tuple, np.ndarray)):
                _write_geom_list(contour, fid, filename)
            else:
                ## Otherwise it is an index and can be written directly to file
                _write_geom(contour, fid, filename)
    
    # If contours is a dictionary, it's just one contour (no loop required)
    else:
        # If it's an array, loop on indexes
        if isinstance(contours['x'], (list, tuple, np.ndarray)):
            _write_geom_list(contours, fid, filename)
        else:
            ## Otherwise it is an index and can be written directly to file
            _write_geom(contours, fid, filename)

    fid.close()

def exp_read(filename):
    """
    Read contours from an exp file.

    This function reads contour data from a file in *.exp format. The function can handle
    files containing multiple contours. Each contour is returned as a dictionary
    containing coordinate data, metadata, and geometric properties.

    Parameters
    ----------
    filename : str
        Path to the input exp file

    Returns
    -------
    contours : list of dict
        List of contour data read from the file. Each contour dictionary contains:
        - 'x' : np.ndarray
            X coordinates of the contour points
        - 'y' : np.ndarray  
            Y coordinates of the contour points
        - 'name' : str
            Name of the contour from the file
        - 'density' : float
            Density value for the contour
        - 'nods' : int
            Number of nodes/points in the contour
        - 'icon' : str, optional
            Icon value from the file, if present
        - 'closed' : bool
            Whether the contour is closed (first and last points are identical)

    Raises
    ------
    IOError
        If the input file does not exist

    Notes
    -----
    The function expects exp files with specific formatting including contour headers,
    point counts, density values, and coordinate data. The function handles variations
    in header spacing (e.g., '# Points Count Value' vs '# Points Count  Value').
    Invalid formatting may cause parsing errors.

    Examples
    --------
    >>> contours = exp_read('input.exp')
    >>> print(f"Read {len(contours)} contours")
    >>> for contour in contours:
    ...     print(f"Contour '{contour['name']}' has {contour['nods']} points")
    """

    # Error checks
    if not os.path.exists(filename):
        raise IOError(f"pyissm.tools.exp.exp_read: File {filename} does not exist.")
    
    # Initialise contours
    contours = []
    contour = None
    
    # Open the file for reading and loop over lines
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip blank lines
            if not line:
                continue
            
            # If Name line, start a new contour
            if line.startswith('## Name:'):
                
                ## Save previous contour if it exists
                if contour is not None:
                    contours.append(contour)
                
                ## Create empty contour
                contour = collections.OrderedDict({
                    'name': line.split('## Name:')[1].strip(),
                    'x': [],
                    'y': [],
                })

            # If Icon line, extract information
            elif line.startswith('## Icon:'):
                contour['icon'] = line.split('## Icon:')[1].strip()

            # If Points Count Value line, read point count and density
            ## NOTE: Some files have '# Points Count Value' and some have '# Points Count  Value'. This handles both.
            elif line.startswith('# Points Count'):
                ## Get next line for point count and density
                nods_line = next(f).strip()

                ## Split the line and extract values
                nods_parts = nods_line.split()
                contour['nods'] = int(nods_parts[0])
                contour['density'] = float(nods_parts[1])

            # If X pos Y pos line, read coordinate data
            elif line.startswith('# X pos Y pos'):
                ## Create empty contour coordinate arrays
                contour['x'] = np.empty(contour['nods'])
                contour['y'] = np.empty(contour['nods'])
                
                ## Read the next 'nods' lines for coordinates
                for i in range(contour['nods']):
                    coord_line = next(f).strip()
                    x_str, y_str = coord_line.split()
                    contour['x'][i] = (float(x_str))
                    contour['y'][i] = (float(y_str))

                ## Check if contour is closed
                contour['closed'] = (
                    contour['nods'] > 1
                    and (contour['x'][-1] == contour['x'][0]) 
                    and (contour['y'][-1] == contour['y'][0])
                )                

        # Append the final contour to the list
        if contour is not None:
            contours.append(contour)

    return contours


def isoline(md, field, value=0.0, output="struct", edges=None, amr=None):
    """
    ISOLINE - construct isovalue lines based on field provided

    Usage:
        contours, edges_tria = isoline(md, field, value=0.0, output="struct", edges=None)

    Supported options (MATLAB parity):
        value  : isoline value (default 0)
        output : "struct" (default), "matrix", "longest", or "filename.exp"
        edges  : if provided, reuse precomputed edges_tria (2D mapping)
        amr    : optional object overriding x/y/elements: amr.MeshX, amr.MeshY, amr.MeshElements

    Returns:
        contours   : list of dicts (exp-struct-like) OR an (N,2) array if output="matrix"
        edges_tria : (nE,3) int array mapping each triangle to its unique-edge ids
    """

    # -----------------------------
    # Process geometry / index arrays
    # -----------------------------
    dim = md.mesh.dimension()

    if dim == 3:
        x = np.asarray(md.mesh.x2d).reshape(-1)
        y = np.asarray(md.mesh.y2d).reshape(-1)
        elements = np.asarray(md.mesh.elements2d)
    else:
        x = np.asarray(md.mesh.x).reshape(-1)
        y = np.asarray(md.mesh.y).reshape(-1)
        elements = np.asarray(md.mesh.elements)

    if amr is not None:
        x = np.asarray(amr.MeshX).reshape(-1)
        y = np.asarray(amr.MeshY).reshape(-1)
        elements = np.asarray(amr.MeshElements)

    # z coordinate (optional)
    if hasattr(md.mesh, "z"):
        z = np.asarray(md.mesh.z).reshape(-1)
    else:
        z = np.zeros_like(x)

    # -----------------------------
    # Checks
    # -----------------------------
    if field is None:
        raise RuntimeError("field provided is empty")

    field = np.asarray(field).reshape(-1)

    if dim == 3:
        if field.size != md.mesh.numberofvertices2d:
            raise RuntimeError("field provided should be of size md.mesh.numberofvertices2d")
    else:
        if field.size != x.size:
            raise RuntimeError("field provided should be of size md.mesh.numberofvertices")

    level = float(value)

    # -----------------------------
    # Convert element indexing to 0-based if needed
    # (MATLAB is 1-based; some pyISSM builds still store 1-based connectivity)
    # -----------------------------
    elements = np.asarray(elements, dtype=int)
    if elements.min() == 1:
        elem = elements - 1
    else:
        elem = elements.copy()

    if elem.ndim != 2 or elem.shape[1] != 3:
        raise RuntimeError("elements array must be (numberofelements, 3)")

    numberofelements = elem.shape[0]
    elementslist = np.arange(numberofelements, dtype=int)

    # -----------------------------
    # Build or reuse unique-edge mapping edges_tria
    # -----------------------------
    if edges is not None:
        edges_tria = np.asarray(edges, dtype=int)
        if edges_tria.shape != (numberofelements, 3):
            raise RuntimeError("provided edges_tria has wrong shape")
    else:
        # edges list: [n1 n2] for each triangle edge
        e12 = elem[:, [0, 1]]
        e23 = elem[:, [1, 2]]
        e31 = elem[:, [2, 0]]
        all_edges = np.vstack([e12, e23, e31])

        # unique on sorted node pairs
        all_edges_sorted = np.sort(all_edges, axis=1)

        # unique rows, and inverse mapping to rebuild per-triangle edge ids
        _, inv = np.unique(all_edges_sorted, axis=0, return_inverse=True)

        # inv has length 3*numberofelements; reshape to (numberofelements,3) in same order
        edges_tria = np.column_stack([
            inv[elementslist],
            inv[elementslist + numberofelements],
            inv[elementslist + 2 * numberofelements],
        ])

    # -----------------------------
    # Segment endpoints per triangle
    # -----------------------------
    Seg1 = elem[:, [0, 1]]
    Seg2 = elem[:, [1, 2]]
    Seg3 = elem[:, [2, 0]]

    Seg1_num = edges_tria[:, 0]
    Seg2_num = edges_tria[:, 1]
    Seg3_num = edges_tria[:, 2]

    # Values at tips (nE,2)
    Data1 = field[Seg1]
    Data2 = field[Seg2]
    Data3 = field[Seg3]

    # Ranges
    Range1 = np.sort(Data1, axis=1)
    Range2 = np.sort(Data2, axis=1)
    Range3 = np.sort(Data3, axis=1)

    # segments that contain level
    pos1 = (Range1[:, 0] < level) & (Range1[:, 1] >= level)
    pos2 = (Range2[:, 0] < level) & (Range2[:, 1] >= level)
    pos3 = (Range3[:, 0] < level) & (Range3[:, 1] >= level)

    poselem12 = pos1 & pos2
    poselem13 = pos1 & pos3
    poselem23 = pos2 & pos3

    poselem = np.where(poselem12 | poselem13 | poselem23)[0]
    numelems = poselem.size

    if numelems == 0:
        # MATLAB: warning(...) and return struct([])
        return [], edges_tria

    # -----------------------------
    # Interpolate segment intersection points
    # -----------------------------
    x1 = np.zeros(numelems)
    x2 = np.zeros(numelems)
    y1 = np.zeros(numelems)
    y2 = np.zeros(numelems)
    z1 = np.zeros(numelems)
    z2 = np.zeros(numelems)

    edge_l = np.zeros((numelems, 2), dtype=int)

    for j, eidx in enumerate(poselem):
        # weights for each segment in element eidx
        # (guard division by zero: if constant along segment, weight becomes nan; those
        #  segments shouldn't be selected by the pos logic unless exactly equals level)
        w1 = (level - Data1[eidx, 0]) / (Data1[eidx, 1] - Data1[eidx, 0])
        w2 = (level - Data2[eidx, 0]) / (Data2[eidx, 1] - Data2[eidx, 0])
        w3 = (level - Data3[eidx, 0]) / (Data3[eidx, 1] - Data3[eidx, 0])

        if poselem12[eidx]:
            a0, a1 = Seg1[eidx, 0], Seg1[eidx, 1]
            b0, b1 = Seg2[eidx, 0], Seg2[eidx, 1]

            x1[j] = x[a0] + w1 * (x[a1] - x[a0])
            y1[j] = y[a0] + w1 * (y[a1] - y[a0])
            z1[j] = z[a0] + w1 * (z[a1] - z[a0])

            x2[j] = x[b0] + w2 * (x[b1] - x[b0])
            y2[j] = y[b0] + w2 * (y[b1] - y[b0])
            z2[j] = z[b0] + w2 * (z[b1] - z[b0])

            edge_l[j, 0] = Seg1_num[eidx]
            edge_l[j, 1] = Seg2_num[eidx]

        elif poselem13[eidx]:
            a0, a1 = Seg1[eidx, 0], Seg1[eidx, 1]
            b0, b1 = Seg3[eidx, 0], Seg3[eidx, 1]

            x1[j] = x[a0] + w1 * (x[a1] - x[a0])
            y1[j] = y[a0] + w1 * (y[a1] - y[a0])
            z1[j] = z[a0] + w1 * (z[a1] - z[a0])

            x2[j] = x[b0] + w3 * (x[b1] - x[b0])
            y2[j] = y[b0] + w3 * (y[b1] - y[b0])
            z2[j] = z[b0] + w3 * (z[b1] - z[b0])

            edge_l[j, 0] = Seg1_num[eidx]
            edge_l[j, 1] = Seg3_num[eidx]

        else:  # poselem23
            a0, a1 = Seg2[eidx, 0], Seg2[eidx, 1]
            b0, b1 = Seg3[eidx, 0], Seg3[eidx, 1]

            x1[j] = x[a0] + w2 * (x[a1] - x[a0])
            y1[j] = y[a0] + w2 * (y[a1] - y[a0])
            z1[j] = z[a0] + w2 * (z[a1] - z[a0])

            x2[j] = x[b0] + w3 * (x[b1] - x[b0])
            y2[j] = y[b0] + w3 * (y[b1] - y[b0])
            z2[j] = z[b0] + w3 * (z[b1] - z[b0])

            edge_l[j, 0] = Seg2_num[eidx]
            edge_l[j, 1] = Seg3_num[eidx]

    # -----------------------------
    # Connect segments into polylines (contours)
    # -----------------------------
    # Convert to python lists for easy deletion like MATLAB code
    edge_l = edge_l.tolist()
    x1 = x1.tolist(); x2 = x2.tolist()
    y1 = y1.tolist(); y2 = y2.tolist()
    z1 = z1.tolist(); z2 = z2.tolist()

    contours = []

    def _find_edge(edge_pairs, edge_id):
        """Return (row, col) if edge_id appears in edge_pairs else (None, None)."""
        for r, (a, b) in enumerate(edge_pairs):
            if a == edge_id:
                return r, 0
            if b == edge_id:
                return r, 1
        return None, None

    while len(edge_l) > 0:
        # start a new contour with the first segment
        e1, e2 = edge_l[0]
        xc = [x1[0], x2[0]]
        yc = [y1[0], y2[0]]
        zc = [z1[0], z2[0]]

        # delete row 0
        del edge_l[0]
        del x1[0]; del x2[0]
        del y1[0]; del y2[0]
        del z1[0]; del z2[0]

        # extend to the "left" by matching e1
        r, c = _find_edge(edge_l, e1)
        while r is not None:
            if c == 0:
                # prepend the other endpoint (segment end 2) to the left
                xc = [x2[r]] + xc
                yc = [y2[r]] + yc
                zc = [z2[r]] + zc
                e1 = edge_l[r][1]
            else:
                xc = [x1[r]] + xc
                yc = [y1[r]] + yc
                zc = [z1[r]] + zc
                e1 = edge_l[r][0]

            # delete that segment
            del edge_l[r]
            del x1[r]; del x2[r]
            del y1[r]; del y2[r]
            del z1[r]; del z2[r]

            r, c = _find_edge(edge_l, e1)

        # extend to the "right" by matching e2
        r, c = _find_edge(edge_l, e2)
        while r is not None:
            if c == 0:
                xc = xc + [x2[r]]
                yc = yc + [y2[r]]
                zc = zc + [z2[r]]
                e2 = edge_l[r][1]
            else:
                xc = xc + [x1[r]]
                yc = yc + [y1[r]]
                zc = zc + [z1[r]]
                e2 = edge_l[r][0]

            del edge_l[r]
            del x1[r]; del x2[r]
            del y1[r]; del y2[r]
            del z1[r]; del z2[r]

            r, c = _find_edge(edge_l, e2)

        # save as exp-struct-like dict
        contour = {
            "x": np.asarray(xc, dtype=float),
            "y": np.asarray(yc, dtype=float),
            "z": np.asarray(zc, dtype=float),
            "name": "",
            "nods": len(xc),
            "density": 1,
            "closed": 0,
        }
        contours.append(contour)

    # -----------------------------
    # Process output formats
    # -----------------------------
    if isinstance(output, str):
        out = output
    else:
        out = "struct"

    if out == "matrix":
        # concatenate x/y with NaN separators
        xs = []
        ys = []
        for c in contours:
            xs.append(c["x"])
            ys.append(c["y"])
            xs.append(np.asarray([np.nan]))
            ys.append(np.asarray([np.nan]))
        if len(xs) == 0:
            return np.zeros((0, 2)), edges_tria
        X = np.concatenate(xs)
        Y = np.concatenate(ys)
        return np.column_stack([X, Y]), edges_tria

    if out == "longest":
        if len(contours) == 0:
            return [], edges_tria
        m = np.argmax([c["nods"] for c in contours])
        return contours[m], edges_tria

    # filename.exp saving (optional)
    if out.endswith(".exp"):
        exp_write(contours, out)
        return contours, edges_tria

    # default: "struct"
    return contours, edges_tria