"""
Functions to execute ISSM commands and manage the execution environment.
"""
from copy import copy
import numpy as np


def marshall(md):
    """
    Marshall the model data for execution.

    Parameters:
    md (ModelData): The model data to be marshalled.

    Returns:
    None
    """

    # If verbose solution is enabled, print the name of the model being marshalled
    # if md.verbose.solution:
    #     print(f'Marshalling for {md.miscellaneous.name}.bin')

    # Open file for binary writing
    try:
        fid = open(md.miscellaneous.name + '.bin', 'wb')
    except IOError as e:
        raise IOError(f"Could not open file {md.miscellaneous.name}.bin for writing: {e}")
    
    # List (and sort) all model classes. Sort simply makes it easier to compare binary files
    model_classes = list(vars(md).keys())
    model_classes.sort()

    # model_classes = ['mask']
    
    for model_class in model_classes:
        # Skip certain classes that do not need marshalling
        if model_class in ['results', 'radaroverlay', 'toolkits', 'cluster', 'private']:
            continue
        
        # Check if the model class has a marshall method
        try:
            callable(getattr(md, model_class).marshall_class)
        except Exception as e:
            # raise TypeError(f'Model class {model_class} does not have a marshall method.')
            print(f"Skipping {model_class} due to error: {e}")
            continue


        # Marshall the model class
        try:
            getattr(md, model_class).marshall_class(prefix = f'md.{model_class}',
                                                    md = md, 
                                                    fid = fid)
        except Exception as e:
            raise RuntimeError(f"Error marshalling model class {model_class}: {e}")
        
    # Write "md.EOF" to make sure that the binary file is not corrupt
    WriteData(fid,
              prefix = 'XXX',
              name = 'md.EOF', 
              data = True, 
              format = 'Boolean')

    #close file
    try:
        fid.close()

    except IOError as e:
        print('marshall error message: could not close \'{}.bin\' file for binary writing due to: {}'.format(md.miscellaneous.name, e))


def WriteData(fid, prefix, *, 
              obj = None,
              fieldname = None,
              data = None,
              name = None,
              format = None,
              mattype = 0,
              timeserieslength = -1,
              scale = None,
              yts = None):
    """WriteData - write model field in binary file using efficient numpy methods"""
    
    # Validate and extract data (same as before)
    if format is None:
        raise ValueError(f"'format' parameter is required")

    if obj is not None:
        if fieldname is None:
            raise ValueError(f"'fieldname' is required when 'obj' is provided")
        if name is None:
            name = f"{prefix}.{fieldname}"
        if data is None:
            data = getattr(obj, fieldname)
    else:
        if data is None:
            raise ValueError(f"Either 'obj'+'fieldname' or 'data' must be provided")
        if name is None:
            raise ValueError(f"'name' is required when using 'data' directly")
    
    # Make copy and apply scaling
    data = copy(data)
    data = _apply_scaling_and_time_series(data, format, timeserieslength, scale, yts)
    
    # Write field identifier
    _write_field_name(fid, name)
    
    # Write data based on format
    _write_data_efficiently(fid, data, format, mattype, name)

def _write_field_name(fid, name):
    """Write field name identifier using numpy"""
    name_bytes = name.encode()
    np.array([len(name_bytes)], dtype=np.int32).tofile(fid)
    fid.write(name_bytes)

def _write_data_efficiently(fid, data, format, mattype, name):
    """Write data using efficient numpy methods"""
    
    if format == 'Boolean':
        _write_boolean_efficient(fid, data, name)
    elif format == 'Integer':
        _write_integer_efficient(fid, data, name)
    elif format == 'Double':
        _write_double_efficient(fid, data, name)
    elif format == 'String':
        _write_string_efficient(fid, data)
    elif format in ['IntMat', 'BooleanMat']:
        _write_int_matrix_efficient(fid, data, format, mattype)
    elif format == 'DoubleMat':
        _write_double_matrix_efficient(fid, data, mattype, name)
    elif format == 'CompressedMat':
        _write_compressed_matrix_efficient(fid, data, mattype, name)
    elif format == 'MatArray':
        _write_matrix_array_efficient(fid, data)
    elif format == 'StringArray':
        _write_string_array_efficient(fid, data)
    else:
        raise TypeError(f'WriteData error: format "{format}" not supported! (field: {name})')

def _write_boolean_efficient(fid, data, name):
    """Write boolean using numpy arrays"""
    try:
        # Write record length and code
        record_info = np.array([4 + 4, FormatToCode('Boolean'), int(data)], dtype=np.int64)
        record_info[:1].astype(np.int64).tofile(fid)  # record length as int64
        record_info[1:].astype(np.int32).tofile(fid)  # code and data as int32
    except Exception as err:
        raise ValueError(f'field {name} cannot be marshaled: {err}')

def _write_integer_efficient(fid, data, name):
    """Write integer using numpy arrays"""
    try:
        # Write record length and code
        np.array([4 + 4], dtype=np.int64).tofile(fid)
        record_data = np.array([FormatToCode('Integer'), int(data)], dtype=np.int32)
        record_data.tofile(fid)
    except Exception as err:
        raise ValueError(f'field {name} cannot be marshaled: {err}')

def _write_double_efficient(fid, data, name):
    """Write double using numpy arrays"""
    try:
        # Write record length and code
        np.array([8 + 4], dtype=np.int64).tofile(fid)
        np.array([FormatToCode('Double')], dtype=np.int32).tofile(fid)
        np.array([float(data)], dtype=np.float64).tofile(fid)
    except Exception as err:
        raise ValueError(f'field {name} cannot be marshaled: {err}')

def _write_string_efficient(fid, data):
    """Write string efficiently"""
    data_bytes = data.encode()
    # Write record length, code, and string length
    np.array([len(data_bytes) + 4 + 4], dtype=np.int64).tofile(fid)
    np.array([FormatToCode('String'), len(data_bytes)], dtype=np.int32).tofile(fid)
    fid.write(data_bytes)

def _preprocess_matrix_efficient(data):
    """Preprocess matrix data ensuring numpy array output"""
    if isinstance(data, (bool, int, float)):
        data = np.array([[data]], dtype=np.float64)
    elif isinstance(data, (list, tuple)):
        data = np.array(data, dtype=np.float64)
        if data.ndim == 1:
            data = data.reshape(-1, 1)
    else:
        data = np.asarray(data, dtype=np.float64)
        if data.ndim == 1:
            if data.size > 0:
                data = data.reshape(-1, 1)
            else:
                data = data.reshape(0, 0)
    
    return data

def _write_double_matrix_efficient(fid, data, mattype, name):
    """Write double matrix using numpy tofile for maximum efficiency"""
    data = _preprocess_matrix_efficient(data)
    
    # Handle NaN matrices
    if data.size == 1 and np.all(np.isnan(data)):
        shape = (0, 0)
        data = np.array([], dtype=np.float64).reshape(0, 0)
    else:
        shape = data.shape
    
    # Calculate record length
    recordlength = 4 + 4 + 8 * np.prod(shape) + 4 + 4
    
    try:
        # Write header information
        header = np.array([recordlength], dtype=np.int64)
        header.tofile(fid)
        
        meta_info = np.array([FormatToCode('DoubleMat'), mattype, shape[0], 
                             shape[1] if len(shape) > 1 else 1], dtype=np.int32)
        meta_info.tofile(fid)
        
        # Write matrix data efficiently - this is the big performance gain!
        if data.size > 0:
            # Ensure C-contiguous for efficient writing
            if not data.flags.c_contiguous:
                data = np.ascontiguousarray(data)
            data.astype(np.float64).tofile(fid)
            
    except Exception as err:
        raise ValueError(f'Field {name} cannot be marshaled: {err}')

def _write_int_matrix_efficient(fid, data, format, mattype):
    """Write integer/boolean matrix efficiently"""
    data = _preprocess_matrix_efficient(data)
    
    # Handle NaN matrices
    if data.size == 1 and np.all(np.isnan(data)):
        shape = (0, 0)
        data = np.array([], dtype=np.float64).reshape(0, 0)
    else:
        shape = data.shape
    
    # Calculate record length
    recordlength = 4 + 4 + 8 * np.prod(shape) + 4 + 4
    
    # Write header
    header = np.array([recordlength], dtype=np.int64)
    header.tofile(fid)
    
    meta_info = np.array([FormatToCode(format), mattype, shape[0], 
                         shape[1] if len(shape) > 1 else 1], dtype=np.int32)
    meta_info.tofile(fid)
    
    # Write data (still as float64 to match original format)
    if data.size > 0:
        if not data.flags.c_contiguous:
            data = np.ascontiguousarray(data)
        data.astype(np.float64).tofile(fid)

def _write_compressed_matrix_efficient(fid, data, mattype, name):
    """Write compressed matrix efficiently"""
    data = _preprocess_matrix_efficient(data)
    
    # Handle NaN matrices
    if data.size == 1 and np.all(np.isnan(data)):
        shape = (0, 0)
        n2 = 0
        data = np.array([], dtype=np.float64).reshape(0, 0)
    else:
        shape = data.shape
        n2 = shape[1] if len(shape) > 1 else 1
    
    # Calculate record length for compressed format
    recordlength = 4 + 4 + 8 + 8 + 1 * (shape[0] - 1) * n2 + 8 * n2 + 4 + 4
    
    try:
        # Write header
        np.array([recordlength], dtype=np.int64).tofile(fid)
        np.array([FormatToCode('CompressedMat'), mattype], dtype=np.int32).tofile(fid)
        
        if shape[0] > 0:
            # Compression logic
            A = data[0:shape[0] - 1]
            offsetA = A.min()
            rangeA = A.max() - offsetA
            
            if rangeA == 0:
                A = A * 0
            else:
                A = (A - offsetA) / rangeA * 255.0
            
            # Write dimensions and compression parameters
            dims_and_params = np.array([shape[0], n2], dtype=np.int32)
            dims_and_params.tofile(fid)
            
            compression_params = np.array([offsetA, rangeA], dtype=np.float64)
            compression_params.tofile(fid)
            
            # Write compressed data
            if A.size > 0:
                A.astype(np.uint8).tofile(fid)
            
            # Write last row as doubles
            if shape[0] > 0:
                last_row = data[shape[0] - 1:shape[0], :]
                last_row.astype(np.float64).tofile(fid)
        else:
            # Empty matrix
            np.array([0, 0], dtype=np.int32).tofile(fid)
            np.array([0.0, 0.0], dtype=np.float64).tofile(fid)
            
    except Exception as err:
        raise ValueError(f'Field {name} cannot be marshaled: {err}')

def _write_matrix_array_efficient(fid, data):
    """Write matrix array efficiently"""
    # Calculate total record length
    recordlength = 4 + 4  # number of records + code
    processed_matrices = []
    
    for matrix in data:
        matrix = _preprocess_matrix_efficient(matrix)
        processed_matrices.append(matrix)
        shape = matrix.shape
        recordlength += 4 * 2 + np.prod(shape) * 8  # dimensions + data
    
    # Write header
    np.array([recordlength], dtype=np.int64).tofile(fid)
    np.array([FormatToCode('MatArray'), len(data)], dtype=np.int32).tofile(fid)
    
    # Write each matrix
    for matrix in processed_matrices:
        shape = matrix.shape
        # Write dimensions
        dims = np.array([shape[0], shape[1] if len(shape) > 1 else 1], dtype=np.int32)
        dims.tofile(fid)
        # Write data
        if matrix.size > 0:
            if not matrix.flags.c_contiguous:
                matrix = np.ascontiguousarray(matrix)
            matrix.astype(np.float64).tofile(fid)

def _write_string_array_efficient(fid, data):
    """Write string array efficiently"""
    # Calculate record length
    recordlength = 4 + 4  # array length + code
    encoded_strings = []
    
    for string in data:
        encoded = string.encode()
        encoded_strings.append(encoded)
        recordlength += 4 + len(encoded)
    
    # Write header
    np.array([recordlength], dtype=np.int64).tofile(fid)
    np.array([FormatToCode('StringArray'), len(data)], dtype=np.int32).tofile(fid)
    
    # Write each string
    for encoded_string in encoded_strings:
        np.array([len(encoded_string)], dtype=np.int32).tofile(fid)
        fid.write(encoded_string)

def _apply_scaling_and_time_series(data, datatype, timeserieslength, scale, yts):
    """Apply scaling and time series transformations (same as before)"""
    if datatype == 'MatArray':
        for i in range(len(data)):
            if scale is not None:
                if np.ndim(data[i]) > 1 and data[i].shape[0] == timeserieslength:
                    data[i][:-1, :] = scale * data[i][:-1, :]
                else:
                    data[i] = scale * data[i]
            if (yts is not None and 
                np.ndim(data[i]) > 1 and data[i].shape[0] == timeserieslength):
                data[i][-1, :] = yts * data[i][-1, :]
    else:
        if scale is not None:
            if np.ndim(data) > 1 and data.shape[0] == timeserieslength:
                data[:-1, :] = scale * data[:-1, :]
            elif isinstance(data, list):
                data = [scale * item for item in data]
            else:
                data = scale * data
        
        if (yts is not None and 
            np.ndim(data) > 1 and data.shape[0] == timeserieslength):
            data[-1, :] = yts * data[-1, :]
    
    return data

def FormatToCode(datatype):
    """Convert format string to integer code"""
    format_codes = {
        'Boolean': 1, 'Integer': 2, 'Double': 3, 'String': 4,
        'BooleanMat': 5, 'IntMat': 6, 'DoubleMat': 7, 
        'MatArray': 8, 'StringArray': 9, 'CompressedMat': 10
    }
    
    if datatype not in format_codes:
        raise IOError(f'FormatToCode error: data type "{datatype}" not supported!')
    
    return format_codes[datatype]

    