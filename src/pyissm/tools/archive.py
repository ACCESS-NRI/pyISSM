"""
This module contains functions to read and write ISSM archive files.
"""

import os
import numpy as np
import collections
import struct

def arch_read(filename, fieldname):
    """
    Read data from an ISSM archive file.

    Parameters
    ----------
    filename : str
        Path to the ISSM archive file.
    fieldname : str
        Name of the field to read from the archive.

    Returns
    -------
    np.ndarray
        Data read from the specified field in the archive.

    Raises
    ------
    FileNotFoundError
        If the archive file does not exist.
    KeyError
        If the specified field is not found in the archive.

    Examples
    --------
    >>> data = arch_read('test101.arch', 'fieldname')
    """

    if not os.path.isfile(filename):
        raise FileNotFoundError(f"Archive file '{filename}' not found.")
    
    with open(filename, 'rb') as fid:
        while True:
            result = _read_field(fid)
            if result is None:
                break

            if result['field_name'] == fieldname:
                return result['data']
    
    raise KeyError(f'Field {fieldname} not found in archive')

def _read_field(fid):
    """
    Read a single field from an ISSM archive file.

    Parameters
    ----------
    fid : file object
        Binary file object opened in read mode.

    Returns
    -------
    dict or None
        Dictionary with keys: field_name, size, data_type, data.
        Returns None on EOF.
    """

    def _read_int(fid):
        """
        Read big-endian 32-bit integer from a binary file.
        """
        return struct.unpack('>i', fid.read(4))[0]
    
    def _read_double(fid):
        """
        Read big-endian 64-bit floating-point value from a binary file.
        """
        return struct.unpack('>d', fid.read(8))[0]

    try:
        # Read the name length & check name
        _read_int(fid)
        check_name = _read_int(fid)
        if check_name != 1:
            raise ValueError('pyissm.tools.archive.arch_read: A string was not present at the start of the archive file')
        
        name_len = _read_int(fid)
        field_name = fid.read(name_len).decode('utf-8')
        
        # Read the data
        _read_int(fid)
        data_type = _read_int(fid)

        if data_type == 2:
            data = _read_double(fid)
            size = '1x1'
            dtype = 'double'
            
        elif data_type == 3:
            rows = _read_int(fid)
            cols = _read_int(fid)
            
            # Read matrix and reshape as float array
            values = np.frombuffer(fid.read(rows * cols * 8), dtype = '>f8')
            data = np.array(values, dtype = float).reshape((cols, rows)).T

            size = 'f{rows}x{cols}'
            dtype = 'vector/matrix'
            
        else:
            raise TypeError(f'Unsupported data type: {data_type}')

        return {
            'field_name': field_name,
            'size': size,
            'data_type': dtype,
            'data': data
        }
    
    except struct.error:
        return None