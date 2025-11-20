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