"""
Tools working with *.exp files for ISSM model domains and contours.
"""

import numpy as np


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

