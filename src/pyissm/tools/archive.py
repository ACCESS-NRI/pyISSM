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

    Parameters:
        filename (str): Path to the ISSM archive file.
        fieldname (str): Name of the field to read from the archive.
    Returns:
        data (np.ndarray): Data read from the specified field in the archive.

    Usage:
        data = arch_read(filename, fieldname)
    """

    if not os.path.isfile(filename):
        raise FileNotFoundError(f"Archive file '{filename}' not found.")
    else:
        fid = open(filename, 'rb')
    
    # Initialize output
    archive_results = []

    # Read first results from file
    result = _read_field(fid)

    while result:
        if fieldname == result['field_name']:
            # Found the data we wanted
            archive_results = result['data']
            break
        # Read next results from file
        result = _read_field(fid)

    # Close file
    fid.close()
    return archive_results

def _read_field(fid):
    """
    Procedure to read a field and return a results list with the following 
    attributes:

        result['field_name']    -> the name of the variable that was just read
        result['size']          -> size (dimensions) of the variable just read
        result['data_type']     -> the type of data that was just read
        result['data']          -> the actual data
    """

    try:
        # first, read the string
        # first read the size and continue reading
        struct.unpack('>i', fid.read(struct.calcsize('>i')))[0]  # name length
        check_name = struct.unpack('>i', fid.read(struct.calcsize('>i')))[0]
        if check_name != 1:
            raise ValueError('archread error : a string was not present at the start of the arch file')
        namelen = struct.unpack('>i', fid.read(struct.calcsize('>i')))[0]
        fieldname = struct.unpack('>{}s'.format(namelen), fid.read(namelen))[0]
        # then, read the data
        # first read the size and continue reading
        struct.unpack('>i', fid.read(struct.calcsize('>i')))[0]  # data length
        data_type = struct.unpack('>i', fid.read(struct.calcsize('>i')))[0]

        if data_type == 2:
            # struct.upack scalar
            data = struct.unpack('>d', fid.read(struct.calcsize('>d')))[0]
        elif data_type == 3:
            rows = struct.unpack('>i', fid.read(struct.calcsize('>i')))[0]
            cols = struct.unpack('>i', fid.read(struct.calcsize('>i')))[0]
            raw_data = np.zeros(shape=(rows, cols), dtype=float)
            for i in range(rows):
                raw_data[i, :] = struct.unpack('>{}d'.format(cols), fid.read(cols * struct.calcsize('>d')))
                # The matrix will be struct.upacked in order and will be filled left -> right by column
                # We need to reshape and transpose the matrix so it can be read correctly
            data = raw_data.reshape(raw_data.shape[::-1]).T
        else:
            raise TypeError("Cannot read data type {}".format(data_type))

        # give additional data to user
        if data_type == 2:
            data_size = '1x1'
            data_type_str = 'double'
        elif data_type == 3:
            data_size = '{0}x{1}'.format(rows, cols)
            data_type_str = 'vector/matrix'

        result = collections.OrderedDict()
        result['field_name'] = fieldname.decode('utf8')
        result['size'] = data_size
        result['data_type'] = data_type_str
        result['data'] = data

    except struct.error as e:
        result = None
        print("result is empty due to", e)

    return result