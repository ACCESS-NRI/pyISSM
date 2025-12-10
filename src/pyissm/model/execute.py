"""
ISSM execution commands to handle model marshalling and execution.
"""

from copy import copy
import numpy as np
import os
import datetime
import warnings
import collections

from pyissm import model, tools


def marshall(md):
    """
    Marshall the model data for execution.

    Parameters
    ----------
    md : object
        The model data object to be marshalled. Must contain miscellaneous.name
        attribute and model classes with marshall_class methods.

    Returns
    -------
    None

    Raises
    ------
    IOError
        If the binary output file cannot be opened or closed.
    RuntimeError
        If an error occurs while marshalling a model class.

    Notes
    -----
    This function iterates through all model classes in the md object and
    calls their marshall_class methods to serialize the data to a binary file.
    Certain classes ('results', 'radaroverlay', 'toolkits', 'cluster', 'private')
    are skipped during marshalling. The function writes an EOF marker at the
    end of the file to ensure data integrity.
    """

    # If verbose solution is enabled, print the name of the model being marshalled
    if md.verbose.solution:
        print(f'Marshalling for {md.miscellaneous.name}.bin')

    # Open file for binary writing
    try:
        fid = open(md.miscellaneous.name + '.bin', 'wb')
    except IOError as e:
        raise IOError(f"Could not open file {md.miscellaneous.name}.bin for writing: {e}")
        
    # Iterate over all model classes and marshall them
    for model_class in md.model_class_names():
        ## Skip certain classes that do not need marshalling
        if model_class in ['results', 'radaroverlay', 'toolkits', 'cluster', 'private']:
            continue

        ## Get the model class object
        obj = getattr(md, model_class)
        
        ## Check if the model class has a marshall method
        try:
            callable(obj.marshall_class)
        except Exception as e:
            print(f"Skipping {model_class} due to error: {e}")
            continue

        ## Marshall the model class
        try:
            obj.marshall_class(fid = fid,
                               prefix = f'md.{model_class}',
                               md = md)
        except Exception as e:
            raise RuntimeError(f"Error marshalling model class {model_class}: {e}")
        
    # Write "md.EOF" to make sure that the binary file is not corrupt
    WriteData(fid,
              prefix = 'XXX',
              name = 'md.EOF', 
              data = True, 
              format = 'Boolean')

    # Close file
    try:
        fid.close()

    except IOError as e:
        print(f'marshall error message: could not close \'{md.miscellaneous.name}.bin\' file for binary writing due to: {e}')


def WriteData(fid, 
              prefix, 
              *, 
              obj = None,
              fieldname = None,
              data = None,
              name = None,
              format = None,
              mattype = 0,
              timeserieslength = -1,
              scale = None,
              yts = None
              ):
    """
    Write model field in binary file.

    Parameters
    ----------
    fid : file object
        File identifier for binary output file.
    prefix : str
        Prefix string for field naming.
    obj : object, optional
        Object containing the field to write. If provided, fieldname must also be specified.
    fieldname : str, optional
        Name of the field attribute in obj. Required when obj is provided.
    data : array_like, optional
        Data to write directly. Required when obj is not provided.
    name : str, optional
        Custom name for the field. If not provided, defaults to "{prefix}.{fieldname}" 
        when using obj, or must be specified when using data directly.
    format : str
        Data format specification for binary writing. Required parameter.
    mattype : int, optional
        Matrix type identifier, by default 0.
    timeserieslength : int, optional
        Length of time series data, by default -1.
    scale : float, optional
        Scaling factor to apply to data before writing, by default None.
    yts : float, optional
        Years to seconds conversion factor, by default None.

    Raises
    ------
    ValueError
        If format is not provided.
        If obj is provided but fieldname is missing.
        If neither obj+fieldname nor data is provided.
        If data is provided directly but name is missing.
        
    Notes
    -----
    The function validates input parameters, extracts data from either an object 
    attribute or direct data input, applies optional scaling and time series 
    transformations, and writes the field to a binary file with appropriate 
    formatting.
    """

    # Validate and extract data
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
    _write_data(fid, data, format, mattype, name)

def _write_field_name(fid, name):
    """
    Write a field name identifier to a binary file.

    This function writes a field name to a binary file by first writing the
    length of the name as a 32-bit integer, followed by the name itself as
    bytes.

    Parameters
    ----------
    fid : file object
        File descriptor or file-like object opened in binary write mode.
    name : str
        The field name to write to the file.

    Returns
    -------
    None

    Notes
    -----
    The function encodes the name as bytes and writes:
    1. The length of the encoded name as a 32-bit integer
    2. The encoded name bytes
    """

    name_bytes = name.encode()
    np.array([len(name_bytes)], dtype=np.int32).tofile(fid)
    fid.write(name_bytes)

def _write_data(fid, data, format, mattype, name):
    """
    Write data to file using efficient numpy methods based on format type.

    This function acts as a dispatcher that routes data writing operations
    to appropriate specialized functions based on the specified format.
    
    Parameters
    ----------
    fid : file object
        File identifier or file object to write data to.
    data : array_like
        The data to be written to the file.
    format : str
        Format specification that determines how data should be written.
        Supported formats include: 'Boolean', 'Integer', 'Double', 'String',
        'IntMat', 'BooleanMat', 'DoubleMat', 'CompressedMat', 'MatArray',
        'StringArray'.
    mattype : str
        Matrix type specification used for matrix formats.
    name : str
        Field name associated with the data being written.

    Raises
    ------
    TypeError
        If the specified format is not supported.

    Notes
    -----
    This function uses numpy methods for efficient data writing operations.
    Different data formats are handled by specialized helper functions.
    """
    
    if format == 'Boolean':
        _write_boolean(fid, data, name)
    elif format == 'Integer':
        _write_integer(fid, data, name)
    elif format == 'Double':
        _write_double(fid, data, name)
    elif format == 'String':
        _write_string(fid, data)
    elif format in ['IntMat', 'BooleanMat']:
        _write_int_matrix(fid, data, format, mattype)
    elif format == 'DoubleMat':
        _write_double_matrix(fid, data, mattype, name)
    elif format == 'CompressedMat':
        _write_compressed_matrix(fid, data, mattype, name)
    elif format == 'MatArray':
        _write_matrix_array(fid, data)
    elif format == 'StringArray':
        _write_string_array(fid, data)
    else:
        raise TypeError(f'WriteData error: format "{format}" not supported! (field: {name})')

def _write_boolean(fid, data, name):
    """
    Write a boolean value to a binary file in ISSM format.
    
    This function writes a boolean value to a binary file in ISSM 
    format with record length, type code, and the boolean data converted
    to integer representation.
    
    Parameters
    ----------
    fid : file object
        File descriptor opened in binary write mode where the data will be written.
    data : bool
        Boolean value to be written to the file.
    name : str
        Name of the field being written, used for error reporting.
        
    Raises
    ------
    ValueError
        If the boolean data cannot be marshaled or written to the file.
        
    Notes
    -----
    The function writes data in the following binary format:
    - Record length (8 bytes, int64)
    - Format code (4 bytes, int32) 
    - Boolean data as integer (4 bytes, int32)
    
    The total record length is 8 bytes (4 bytes for code + 4 bytes for data).
    """

    try:
        # Write record length and code
        record_info = np.array([4 + 4, format_to_code('Boolean'), int(data)], dtype=np.int64)
        record_info[:1].astype(np.int64).tofile(fid)  # record length as int64
        record_info[1:].astype(np.int32).tofile(fid)  # code and data as int32
    except Exception as err:
        raise ValueError(f'field {name} cannot be marshaled: {err}')

def _write_integer(fid, data, name):
    """
    Write an integer value to a binary file in ISSM format.

    This function writes an integer value to a binary file using a structured format
    that includes record length, type code, and the actual integer data.

    Parameters
    ----------
    fid : file object
        File identifier/handle opened in binary write mode where the integer
        data will be written.
    data : int or int-like
        The integer value to be written to the file. Will be converted to int
        if not already an integer type.
    name : str
        Name or identifier of the field being written, used for error reporting
        and debugging purposes.

    Raises
    ------
    ValueError
        If the field cannot be marshaled due to data conversion issues or
        file writing errors. The original exception is included in the error message.

    Notes
    -----
    The function writes data in the following binary format:
    1. Record length (8 bytes, int64): Total size of the following data (8 bytes)
    2. Type code (4 bytes, int32): Format code for 'Integer' type
    3. Integer value (4 bytes, int32): The actual integer data
    """

    try:
        # Write record length and code
        np.array([4 + 4], dtype=np.int64).tofile(fid)
        record_data = np.array([format_to_code('Integer'), int(data)], dtype=np.int32)
        record_data.tofile(fid)
    except Exception as err:
        raise ValueError(f'field {name} cannot be marshaled: {err}')

def _write_double(fid, data, name):
    """
    Write a double precision floating point value to a binary file.

    This function writes a double precision value to a binary file in ISSM
    format that includes a record length, format code, and the actual data value.

    Parameters
    ----------
    fid : file object
        File handle opened in binary write mode where the data will be written.
    data : float or numeric
        The numeric value to be written as a double precision float.
    name : str
        Name identifier for the field being written, used for error reporting.

    Raises
    ------
    ValueError
        If the data cannot be marshaled or written to the file, with details
        about which field caused the error.

    Notes
    -----
    The function writes data in the following binary format:
    1. Record length (8 + 4 bytes) as int64
    2. Format code for 'Double' type as int32  
    3. The actual data value as float64
    """

    try:
        # Write record length and code
        np.array([8 + 4], dtype=np.int64).tofile(fid)
        np.array([format_to_code('Double')], dtype=np.int32).tofile(fid)
        np.array([float(data)], dtype=np.float64).tofile(fid)
    except Exception as err:
        raise ValueError(f'field {name} cannot be marshaled: {err}')

def _write_string(fid, data):
    """
    Write string data to a binary file in a structured format.

    This function writes string data to a binary file in ISSM format that includes
    record length, format code, string length, and the actual string data as bytes.

    Parameters
    ----------
    fid : file-like object
        File identifier or file handle opened in binary write mode.
    data : str
        The string data to be written to the file.

    Notes
    -----
    The function writes data in the following binary format:
    1. Record length (8 bytes, int64): Total length of code + string length + string data
    2. Format code (4 bytes, int32): Code identifying this as a string type
    3. String length (4 bytes, int32): Length of the encoded string in bytes
    4. String data (variable length): UTF-8 encoded string bytes

    The record length is calculated as: len(encoded_string) + 4 + 4 bytes.
    """

    data_bytes = data.encode()
    # Write record length, code, and string length
    np.array([len(data_bytes) + 4 + 4], dtype=np.int64).tofile(fid)
    np.array([format_to_code('String'), len(data_bytes)], dtype=np.int32).tofile(fid)
    fid.write(data_bytes)

def _preprocess_matrix(data):
    """
    Preprocess matrix data ensuring numpy array output.

    This function converts various input data types into a standardized numpy array
    format suitable for matrix operations. Scalar values are converted to 2D arrays,
    1D arrays are reshaped to column vectors, and existing arrays are converted to
    float64 dtype.
    
    Parameters
    ----------
    data : bool, int, float, list, tuple, or array-like
        Input data to be preprocessed. Can be:
        - Scalar values (bool, int, float): converted to 2D array
        - Sequences (list, tuple): converted to numpy array
        - Array-like objects: converted to numpy array with float64 dtype

    Returns
    -------
    numpy.ndarray
        Preprocessed data as a numpy array with dtype float64.
        - Scalars become 1x1 arrays
        - 1D arrays with size > 0 become column vectors (n, 1)
        - Empty 1D arrays become (0, 0) arrays
        - Higher dimensional arrays preserve their shape

    Notes
    -----
    All output arrays have dtype numpy.float64 regardless of input type.
    One-dimensional arrays are always reshaped to ensure proper matrix dimensions.
    """

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

def _write_double_matrix(fid, data, mattype, name):
    """
    Write double matrix to binary file.

    This function writes a double-precision matrix to a binary file in a specific
    format that includes header information and matrix data. It handles special
    cases like NaN matrices and ensures optimal performance through efficient
    memory layout.

    Parameters
    ----------
    fid : file-like object
        File handle opened in binary write mode where the matrix will be written.
    data : array_like
        Input matrix data to be written. Will be converted to numpy array and
        cast to float64 precision.
    mattype : int
        Matrix type identifier used in the binary format header.
    name : str
        Name of the matrix field, used for error reporting purposes.

    Raises
    ------
    ValueError
        If the matrix data cannot be marshaled to the binary format, with
        details about the specific error encountered.

    Notes
    -----
    The binary format written includes:
    - Record length (8 bytes, int64)
    - Format code for 'DoubleMat' (4 bytes, int32)  
    - Matrix type (4 bytes, int32)
    - Number of rows (4 bytes, int32)
    - Number of columns (4 bytes, int32)
    - Matrix data in row-major order (8*rows*cols bytes, float64)
    Special handling is performed for NaN matrices, which are written as
    empty (0x0) matrices. Data is automatically converted to C-contiguous
    layout for optimal write performance.
    """

    data = _preprocess_matrix(data)

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
        
        meta_info = np.array([format_to_code('DoubleMat'), mattype, shape[0], 
                             shape[1] if len(shape) > 1 else 1], dtype=np.int32)
        meta_info.tofile(fid)
        
        # Write matrix data
        if data.size > 0:
            # Ensure C-contiguous for efficient writing
            if not data.flags.c_contiguous:
                data = np.ascontiguousarray(data)
            data.astype(np.float64).tofile(fid)
            
    except Exception as err:
        raise ValueError(f'Field {name} cannot be marshaled: {err}')

def _write_int_matrix(fid, data, format, mattype):
    """
    Write integer/boolean matrix to a binary file.

    This function handles the writing of integer or boolean matrices
    to a binary file with a specific format including headers and metadata.
    Matrices containing NaN values are handled as empty matrices.

    Parameters
    ----------
    fid : file object
        File descriptor or file-like object opened in binary write mode.
    data : array_like
        Input matrix data to be written. Can be integer or boolean type.
        NaN matrices are converted to empty arrays.
    format : str or int
        Format specification that will be converted to a format code.
    mattype : int
        Matrix type identifier to be written in the metadata.

    Notes
    -----
    The function writes data in the following binary format:
    1. Record length (int64)
    2. Format code (int32)
    3. Matrix type (int32) 
    4. Number of rows (int32)
    5. Number of columns (int32)
    6. Matrix data as float64 values
    Data is converted to C-contiguous layout if necessary and cast to float64
    for consistency with the original format specification.
    For matrices containing only NaN values, an empty (0,0) matrix is written
    instead of the NaN data.
    """

    data = _preprocess_matrix(data)

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
    
    meta_info = np.array([format_to_code(format), mattype, shape[0], 
                         shape[1] if len(shape) > 1 else 1], dtype=np.int32)
    meta_info.tofile(fid)
    
    # Write data
    if data.size > 0:
        if not data.flags.c_contiguous:
            data = np.ascontiguousarray(data)
        data.astype(np.float64).tofile(fid)

def _write_compressed_matrix(fid, data, mattype, name):
    """Write a matrix to file using compressed format.
    This function writes a matrix to a binary file using a custom compressed format.
    The compression scheme stores all rows except the last as compressed uint8 values,
    while the last row is stored as uncompressed float64 values. NaN matrices are
    handled as special empty cases.

    Parameters
    ----------
    fid : file object
        File descriptor opened in binary write mode where the matrix will be written.
    data : array_like
        Input matrix data to be compressed and written. Can be any array-like structure
        that can be converted to a numpy array.
    mattype : int
        Matrix type identifier used in the file format header.
    name : str
        Name of the matrix field, used for error reporting.

    Raises
    ------
    ValueError
        If the matrix cannot be marshaled or written to file, with details about
        the specific field that failed.

    Notes
    -----
    The compressed format uses the following structure:
    - Record length (int64)
    - Format code and matrix type (2x int32)
    - Matrix dimensions (2x int32)
    - Compression parameters: offset and range (2x float64)
    - Compressed matrix data excluding last row (uint8)
    - Last row as uncompressed doubles (float64)
    The compression algorithm normalizes values to 0-255 range using:
    compressed = (original - offset) / range * 255
    Special handling is provided for NaN matrices, which are written as empty
    matrices with zero dimensions.
    """

    data = _preprocess_matrix(data)

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
        np.array([format_to_code('CompressedMat'), mattype], dtype=np.int32).tofile(fid)
        
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

def _write_matrix_array(fid, data):
    """
    Write an array of matrices to a binary file in ISSM format.

    This function efficiently writes multiple matrices to a binary file using
    the ISSM binary format. Each matrix is preprocessed and written with its
    dimensions and data in a contiguous format.

    Parameters
    ----------
    fid : file object
        Binary file object opened for writing where the matrix array will be written.
    data : array_like
        Sequence of matrices to write. Each matrix will be preprocessed and converted
        to float64 format before writing.

    Notes
    -----
    The binary format includes:
    - Header with total record length (int64)
    - Format code and number of matrices (int32)
    - For each matrix: dimensions (2 int32 values) followed by data (float64)
    All matrices are converted to C-contiguous arrays and float64 dtype before
    writing to ensure consistent binary format.

    See Also
    --------
    _preprocess_matrix : Function used to preprocess each matrix before writing
    format_to_code : Function to convert format string to numeric code
    """

    # Calculate total record length
    recordlength = 4 + 4  # number of records + code
    processed_matrices = []
    
    for matrix in data:
        matrix = _preprocess_matrix(matrix)
        processed_matrices.append(matrix)
        shape = matrix.shape
        recordlength += 4 * 2 + np.prod(shape) * 8  # dimensions + data
    
    # Write header
    np.array([recordlength], dtype=np.int64).tofile(fid)
    np.array([format_to_code('MatArray'), len(data)], dtype=np.int32).tofile(fid)
    
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

def _write_string_array(fid, data):
    """
    Write string array to binary file in ISSM format.

    This function writes an array of strings to a binary file using the ISSM
    binary format. Each string is encoded and written with its length prefix.

    Parameters
    ----------
    fid : file-like object
        File descriptor or file-like object opened in binary write mode.
    data : array-like of str
        Array or list of strings to write to the file.

    Notes
    -----
    The binary format consists of:
    - Record length (int64): Total size of the record in bytes
    - Format code (int32): Identifier for StringArray format
    - Array length (int32): Number of strings in the array
    - For each string:
        - String length (int32): Number of bytes in the encoded string
        - String data (bytes): UTF-8 encoded string content
    The record length includes the format code, array length, all string
    length prefixes, and the encoded string data.
    """

    # Calculate record length
    recordlength = 4 + 4  # array length + code
    encoded_strings = []
    
    for string in data:
        encoded = string.encode()
        encoded_strings.append(encoded)
        recordlength += 4 + len(encoded)
    
    # Write header
    np.array([recordlength], dtype=np.int64).tofile(fid)
    np.array([format_to_code('StringArray'), len(data)], dtype=np.int32).tofile(fid)
    
    # Write each string
    for encoded_string in encoded_strings:
        np.array([len(encoded_string)], dtype=np.int32).tofile(fid)
        fid.write(encoded_string)

def _apply_scaling_and_time_series(data, datatype, timeserieslength, scale, yts):
    """
    Apply scaling and time series transformations to data arrays.

    This function applies scaling factors and year-to-second (yts) conversions
    to input data, handling both single arrays and MatArray collections. For
    time series data, scaling is applied to all but the last row, while yts
    conversion is applied only to the last row (where the time variable is located).

    Parameters
    ----------
    data : array_like or list
        Input data to be scaled. Can be a numpy array, list, or collection
        of arrays when datatype is 'MatArray'.
    datatype : str
        Type of data structure. If 'MatArray', data is treated as a collection
        of arrays; otherwise as a single data structure.
    timeserieslength : int
        Expected length of time series dimension. Used to identify time series
        data for special handling.
    scale : float or None
        Scaling factor to apply to data. If None, no scaling is performed.
        For time series data, applied to all rows except the last.
    yts : float or None
        Year-to-second conversion factor. If None, no conversion is performed.
        For time series data, applied only to the last row (where the time variable is located).

    Returns
    -------
    data : array_like or list
        Transformed data with scaling and time series conversions applied.
        Type matches the input data type.

    Notes
    -----
    - For MatArray datatype, each element in the collection is processed
      independently
    - Time series data is identified by having more than 1 dimension and
      first dimension equal to timeserieslength
    - Scaling is applied element-wise using numpy broadcasting rules
    - The function modifies data in-place for arrays but returns new objects
      for lists
    """

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

def format_to_code(datatype):
    """
    Convert format string to integer code.

    This function maps data type strings to their corresponding integer codes
    used for data serialization and format identification.

    Parameters
    ----------
    datatype : str
        The data type string to convert. Supported types are:
        'Boolean', 'Integer', 'Double', 'String', 'BooleanMat', 
        'IntMat', 'DoubleMat', 'MatArray', 'StringArray', 'CompressedMat'.

    Returns
    -------
    int
        The integer code corresponding to the input data type.

    Raises
    ------
    IOError
        If the provided datatype is not supported.

    Examples
    --------
    >>> format_to_code('Double')
    3
    >>> format_to_code('StringArray')
    9
    """

    format_codes = {
        'Boolean': 1,
        'Integer': 2,
        'Double': 3,
        'String': 4,
        'BooleanMat': 5,
        'IntMat': 6,
        'DoubleMat': 7,
        'MatArray': 8,
        'StringArray': 9,
        'CompressedMat': 10
    }
    
    if datatype not in format_codes:
        raise IOError(f'format_to_code error: data type "{datatype}" not supported!')
    
    return format_codes[datatype]

def solve(md,
          solution_string,
          batch = False,
          check_consistency = True,
          restart = None,
          load_only = False,
          no_log = False,
          runtime_name = True):
    """
    Solve an ISSM model with the specified solution type.

    This function configures and executes a solution for an ISSM Model.
    It maps solution strings to their corresponding solution classes, performs model
    consistency checks, and manages runtime naming conventions.

    Parameters
    ----------
    md : object
        Model data structure containing the ISSM model configuration.
    solution_string : str
        String identifier for the solution type. Supported values include:
        - 'sb', 'stressbalance' : Stress balance solution
        - 'mt', 'masstransport' : Mass transport solution  
        - 'oceant', 'oceantransport' : Ocean transport solution
        - 'th', 'thermal' : Thermal solution
        - 'st', 'steadystate' : Steady state solution
        - 'tr', 'transient' : Transient solution
        - 'mc', 'balancethickness' : Balance thickness solution
        - 'mcsoft' : Balance thickness soft solution
        - 'bv', 'balancevelocity' : Balance velocity solution
        - 'bsl', 'bedslope' : Bed slope solution
        - 'ssl', 'surfaceslope' : Surface slope solution
        - 'hy', 'hydrology' : Hydrology solution
        - 'da', 'damageevolution' : Damage evolution solution
        - 'gia' : GIA solution
        - 'lv', 'love' : Love solution
        - 'esa' : ESA solution
        - 'smp', 'sampling' : Sampling solution
    batch : bool, optional
        Whether to run in batch mode, by default False.
    check_consistency : bool, optional
        Whether to perform model consistency checks before solving, by default True.
    restart : str, optional
        Directory name (relative to execution directory) where restart file is located, by default None.
    load_only : bool, optional
        Whether to only load the solution without executing, by default False.
    no_log : bool, optional
        Whether to suppress loading logs, by default False.
    runtime_name : bool, optional
        Whether to generate a unique runtime name based on timestamp and process ID,
        by default True. If False, uses md.miscellaneous.name as runtime name.

    Raises
    ------
    ValueError
        If the provided solution_string is not recognized or supported.

    Notes
    -----
    The function modifies the model structure in-place by setting the solution type
    and runtime name. When runtime_name is True, a unique identifier is generated
    using the current timestamp and process ID to avoid clobbering existing runs.

    Examples
    --------
    >>> solve(md, 'stressbalance')
    >>> solve(md, 'transient', check_consistency = False)
    >>> solve(md, 'thermal', runtime_name = False)
    """

    # Map solution string to solution names
    solution_map = {
        'sb': 'StressbalanceSolution', 'stressbalance': 'StressbalanceSolution',
        'mt': 'MasstransportSolution', 'masstransport': 'MasstransportSolution',
        'oceant': 'OceantransportSolution', 'oceantransport': 'OceantransportSolution',
        'th': 'ThermalSolution', 'thermal': 'ThermalSolution',
        'st': 'SteadystateSolution', 'steadystate': 'SteadystateSolution',
        'tr': 'TransientSolution', 'transient': 'TransientSolution',
        'mc': 'BalancethicknessSolution', 'balancethickness': 'BalancethicknessSolution',
        'mcsoft': 'BalancethicknessSoftSolution',
        'bv': 'BalancevelocitySolution', 'balancevelocity': 'BalancevelocitySolution',
        'bsl': 'BedSlopeSolution', 'bedslope': 'BedSlopeSolution',
        'ssl': 'SurfaceSlopeSolution', 'surfaceslope': 'SurfaceSlopeSolution',
        'hy': 'HydrologySolution', 'hydrology': 'HydrologySolution',
        'da': 'DamageEvolutionSolution', 'damageevolution': 'DamageEvolutionSolution',
        'gia': 'GiaSolution',
        'lv': 'LoveSolution', 'love': 'LoveSolution',
        'esa': 'EsaSolution',
        'smp': 'SamplingSolution', 'sampling': 'SamplingSolution',
    }

    ## Check that solution string is valid and extract solution name
    key = solution_string.lower()
    if key not in solution_map:
        raise ValueError(f'solve error: solution "{solution_string}" not recognized!')
    solution = solution_map[key]
    
    # Parse options and fields
    ## Set solution in model structure
    md.private.solution = solution

    ## Check model consistency
    if check_consistency:
        if md.verbose.solution:
            print('Checking model consistency...')
        is_model_self_consistent(md)
    
    ## If using restart, use the provided runtime name
    if restart is not None:
        md.private.runtimename = restart

    ## If runtime_name is true, generate a unique runtime name
    if runtime_name:
        now = datetime.datetime.now()
        md.private.runtimename = (
            f"{md.miscellaneous.name}-{now.month:02d}-{now.day:02d}-"
            f"{now.year:04d}-{now.hour:02d}-{now.minute:02d}-"
            f"{now.second:02d}-{os.getpid()}"
        )
    else:
        ### Otherwise use the model name as runtime name
        md.private.runtimename = md.miscellaneous.name

    ## Error when using md.settings.io_gather = 0
    ## NOTE: parse_results_from_disk_io_split is not yet implemented in pyISSM.
    if md.settings.io_gather == 0:
        raise NotImplementedError('pyissm.execute.solve: Reading model results when using md.settings.io_gather = 0 is not yet supported.')

    ## If running QMU analysis, some preprocessing of Dakota files is needed
    if md.qmu.isdakota:
        if md.verbose.solution:
            print('Preprocessing dakota files...')
        preprocess_qmu(md)

    ## If load_only is true, skip the actual solve
    if load_only:
        if md.verbose.solution:
            print('Loading results from cluster...')
        md = load_results_from_cluster(md, no_log = no_log)
        return md
    
    # Write all input files (.bin, .toolkits, build queue script)
    ## Extract model_name
    model_name = md.miscellaneous.name

    ## Marshall model (write .bin file)
    marshall(md)

    ## Write toolkits file
    md.toolkits.write_toolkits_file(model_name + '.toolkits')

    ## Build queue script (and associated logs, if necessary)
    md.cluster.build_queue_script(
        dir_name = md.private.runtimename,
        model_name = model_name,
        solution = md.private.solution,
        io_gather = md.settings.io_gather,
        is_valgrind = md.debug.valgrind,
        is_gprof = md.debug.gprof,
        is_dakota = md.qmu.isdakota,
        is_ocean_coupling = md.transient.isoceancoupling
    )
    
    # Upload all required files (if no restart is provided)
    if restart is None:
        ## Create list of files to upload
        file_list = [model_name + ext for ext in ['.bin', '.toolkits']]
    
        ## Append appropriate queue script to file_list
        if tools.config.is_pc():
            file_list.append(model_name + '.bat')
        else:
            file_list.append(model_name + '.queue')

        ## If using dakota, append dakota input file to file_list
        if md.qmu.isdakota:
            file_list.append('qmu.in')

        ## Upload all files to cluster
        print('Uploading files to cluster...')
        md.cluster.upload_queue_job(model_name, md.private.runtimename, file_list)

    # Launch job
    md.cluster.launch_queue_job(
        model_name,
        md.private.runtimename,
        restart,
        batch
        )
    
    # Early return if using batch
    if batch:
        if md.verbose.solution:
            print('batch mode: not launching job interactively')
            print('launch solution sequence on remote cluster by hand')
        return md
    
    # Wait for job to finish
    if md.settings.waitonlock > 0:
        if md.verbose.solution:
            print('Waiting for job to complete...')
        done = wait_on_lock(md)
        if md.verbose.solution:
            print('Job completed -- loading results from cluster...')
            md = load_results_from_cluster(md)
    else:
        print('Model results must be loaded manually with md = load_results_from_cluster(md)')

    return md


def _get_analysis_for_solution(solution_type):
    """
    Map solution type to analysis type.

    This function returns the analysis type corresponding to a given solution type.
    It is used to determine the appropriate analysis method for different solution
    strategies in the ISSM framework.

    Parameters
    ----------
    solution_type : str
        The solution type string. Supported values include:


    Returns
    -------
    str
        The corresponding analysis type string.

    Raises
    ------
    ValueError
        If the provided solution_type is not recognized.
    """

    analyses_map = {'StressbalanceSolution': ['StressbalanceAnalysis',
                                            'StressbalanceVerticalAnalysis',
                                            'StressbalanceSIAAnalysis',
                                            'L2ProjectionBaseAnalysis'],
                    'SteadystateSolution': ['StressbalanceAnalysis',
                                            'StressbalanceVerticalAnalysis',
                                            'StressbalanceSIAAnalysis',
                                            'L2ProjectionBaseAnalysis',
                                            'ThermalAnalysis',
                                            'MeltingAnalysis',
                                            'EnthalpyAnalysis',
                                            'AgeAnalysis'],
                    'ThermalSolution': ['EnthalpyAnalysis',
                                        'ThermalAnalysis',
                                        'MeltingAnalysis'],
                    'MasstransportSolution': ['MasstransportAnalysis'],
                    'OceantransportSolution': ['OceantransportAnalysis'],
                    'BalancethicknessSolution': ['BalancethicknessAnalysis'],
                    'Balancethickness2Solution': ['Balancethickness2Analysis'],
                    'BalancethicknessSoftSolution': ['BalancethicknessAnalysis'],
                    'BalancevelocitySolution': ['BalancevelocityAnalysis'],
                    'SurfaceSlopeSolution': ['L2ProjectionBaseAnalysis'],
                    'BedSlopeSolution': ['L2ProjectionBaseAnalysis'],
                    'GiaSolution': ['GiaIvinsAnalysis'],
                    'LoveSolution': ['LoveAnalysis'],
                    'EsaSolution': ['EsaAnalysis'],
                    'TransientSolution': ['StressbalanceAnalysis',
                                        'StressbalanceVerticalAnalysis',
                                        'StressbalanceSIAAnalysis',
                                        'L2ProjectionBaseAnalysis',
                                        'ThermalAnalysis',
                                        'MeltingAnalysis',
                                        'EnthalpyAnalysis',
                                        'MasstransportAnalysis',
                                        'OceantransportAnalysis',
                                        'HydrologyShaktiAnalysis',
                                        'HydrologyGladsAnalysis',
                                        'HydrologyShreveAnalysis',
                                        'HydrologyTwsAnalysis',
                                        'HydrologyDCInefficientAnalysis',
                                        'HydrologyDCEfficientAnalysis',
                                        'SealevelchangeAnalysis',
                                        'AgeAnalysis',
                                        'HydrologyArmapwAnalysis',
                                        'DebrisAnalysis'],
                        'SealevelchangeSolution': ['SealevelchangeAnalysis'],
                        'HydrologySolution': ['L2ProjectionBaseAnalysis',
                                            'HydrologyShreveAnalysis',
                                            'HydrologyDCInefficientAnalysis',
                                            'HydrologyDCEfficientAnalysis',
                                            'HydrologyGladsAnalysis',
                                            'HydrologyShaktiAnalysis',
                                            'HydrologyTwsAnalysis',
                                            'HydrologyArmapwAnalysis'],
                        'DamageEvolutionSolution': ['DamageEvolutionAnalysis'],
                        'SamplingSolution': ['SamplingAnalysis']
                        }

    try:
        return analyses_map[solution_type]
    except KeyError:
        raise TypeError(f"pyissm.execute._get_analysis_for_solution: Solution type '{solution_type}' not supported.")


def is_model_self_consistent(md):
    """
    Check that all relevant model classes in the given model are self-consistent.
    This function iterates over the model classes returned by md.model_class_names(),
    skips a small set of classes that are not subject to consistency checking
    ('results', 'debug', 'radaroverlay'), and invokes a consistency check on each
    model class that provides a check_consistency method. The routine:
    - Initializes md.private.isconsistent to True.
    - Retrieves the model solution and any associated analyses.
    - For each model class:
        - If the class does not expose a callable `check_consistency` it is skipped
            (a message is printed to stdout).
        - Otherwise, calls obj.check_consistency(md, solution, analyses).
        - If that call raises an exception, records a message via md.check_message
            and sets md.private.isconsistent to False.
    - If any consistency check fails (md.private.isconsistent becomes False),
        raises a RuntimeError to indicate the overall model is not consistent.

    Parameters
    ----------
    md : object
            The model object to validate. The object is expected to provide:
            - md.private, an attribute with at least `isconsistent` (bool) and
                `solution` attributes;
            - md.model_class_names(), a callable that returns an iterable of model
                class name strings;
            - for each name returned by md.model_class_names(), an attribute on `md`
                with that name, representing the model class object;
            - md.check_message(msg: str), a callable used to record consistency
                failure messages.

    Returns
    -------
    None
            This function performs in-place updates to `md` (notably
            md.private.isconsistent) and does not return a value.

    Raises
    ------
    RuntimeError
            If any per-class consistency check fails. In that case, md.private.isconsistent
            will be set to False and a RuntimeError is raised after the failing check(s)
            so callers are made aware that the model is not consistent.

    Notes
    -----
    - Classes named 'results', 'debug', and 'radaroverlay' are intentionally skipped.
    - If a model class does not provide a callable `check_consistency` attribute,
        a message will be printed to stdout and the class will be skipped.
    - Consistency checks are expected to follow the signature:
        check_consistency(md, solution, analyses). Any exception raised by a
        check_consistency implementation will be caught, recorded via
        md.check_message, and treated as a failure.
    - This function mutates md.private.isconsistent; callers can inspect that value
        to determine consistency state prior to handling the RuntimeError.
    """

    # Initialize consistency flag
    md.private.isconsistent = True

    # Get solution and associated analyses
    solution = md.private.solution
    analyses = _get_analysis_for_solution(solution)

    for model_class in md.model_class_names():
        # Skip select model classes
        if model_class in {'results', 'debug', 'radaroverlay'}:
            continue

        # Get model class object
        obj = getattr(md, model_class)

        # Check if the model class has a check_consistency method
        try:
            callable(obj.check_consistency)
        except Exception as e:
            # NOTE: using md.check_message here would throw an inconsistency error since it modifies md.private.isconsistent. Just skip these classes for now.
            print(f"No check_consistency() available for model class '{model_class}'.")
            continue

        # Perform consistency check
        try:
            obj.check_consistency(md, solution, analyses)
        except Exception as e:
            md.check_message(f"Consistency check failed for model class '{model_class}': {e}")
            md.private.isconsistent = False

        # If an inconsistency is found, md.private.isconsistent is set to False
        if not md.private.isconsistent:
            raise RuntimeError("pyissm.execute.is_model_self_consistent: model not consistent — see messages above.")
    
    return

def preprocess_qmu(md):
    print('preprocess_qmu not implemented yet')
    return

def postprocess_qmu(md):
    print('postprocess_qmu not implemented yet')
    return

def wait_on_lock(md):
    print('wait_on_lock not implemented yet')
    return

def load_results_from_cluster(md,
                              no_log = False,
                              runtime_name = None):
    """
    Load results from cluster after job completion.

    Downloads output files from a remote cluster, loads the results into the model,
    and cleans up temporary files. Handles both standard ISSM runs and Dakota
    uncertainty quantification runs.

    Parameters
    ----------
    md : model
        ISSM model object containing cluster configuration and run parameters.
    no_log : bool, optional
        If True, skip downloading log files (.outlog and .errlog).
        Default is False.
    runtime_name : str, optional
        Name of the runtime directory on the cluster. If None, uses
        md.private.runtimename. Default is None.

    Returns
    -------
    model
        Updated model object with loaded results from the cluster run.

    Notes
    -----
    This function performs the following operations:
    1. Downloads output files from the cluster including:
    - Log files (.outlog, .errlog) unless no_log=True
    - Binary output files (.outbin or multiple files for Dakota sampling)
    - Dakota-specific files (qmu.err, qmu.out, dakota_tabular.dat, .stats)
    2. Loads results from the binary output file into the model object
    3. Cleans up downloaded files and removes temporary files
    4. If run on the same platform as the cluster, removes input files
    (.bin, .toolkits, qmu.in, .queue/.bat)
    The function handles different file patterns depending on whether Dakota
    uncertainty quantification is used and what type of Dakota method is employed.
    Warnings are issued if expected files are not found during cleanup.

    """

    # Set default runtime_name
    if runtime_name is None:
        runtime_name = md.private.runtimename
    
    # Create list of files to download
    if no_log:
        file_list = []
    else:
        file_list = [md.miscellaneous.name + '.outlog', md.miscellaneous.name + '.errlog']

    # If using dakota, append dakota output files to file_list
    if md.qmu.isdakota:
        file_list.append('qmu.err')
        file_list.append('qmu.out')
        if 'tabular_graphics_data' in list(vars(md.qmu).keys()):
            file_list.append('dakota_tabular.dat')
        if md.qmu.output and md.qmu.statistics.method[0]['name'] == 'None':
            if md.qmu.method.method == 'nond_sampling':
                for i in range(md.qmu.method.params.samples):
                    file_list.append(md.miscellaneous.name + '.outbin.' + str(i + 1))
        if md.qmu.statistics.method[0]['name'] != 'None':
            file_list.append(md.miscellaneous.name + '.stats')
    else:
        ## Otherwise, append standard output binary file
        file_list.append(md.miscellaneous.name + '.outbin')
    
    # Download all files from cluster
    md.cluster.download(md.private.runtimename, file_list)

    # Load results from downloaded files
    md = load_results_from_disk(md, md.miscellaneous.name + '.outbin')

    # Erase log and output files
    for i in range(len(file_list)):
        try:
            os.remove(file_list[i])
        except OSError:
            warnings.warn(f'pyissm.execute.load_results_from_cluster: File {file_list[i]} not found.')
    
    # Remove initial tar.gz file
    if not tools.config.is_pc() and os.path.exists(md.private.runtimename + '.tar.gz'):
        os.remove(md.private.runtimename + '.tar.gz')
    
    # Erase all input files if run was carried out on same platform
    hostname = tools.config.get_hostname()
    if hostname.lower() == md.cluster.name.lower():
        
        ## Remove bin and toolkits files
        os.remove(md.miscellaneous.name + '.bin')
        os.remove(md.miscellaneous.name + '.toolkits')
        
        ## If using dakota, remove dakota input file
        if md.qmu.isdakota and os.path.exists('qmu.in'):
            os.remove('.qmu.in')

        ## Remove queue script
        if not tools.config.is_pc():
            if os.path.exists(md.miscellaneous.name + '.queue'):
                os.remove(md.miscellaneous.name + '.queue')
        else:
            if os.path.exists(md.miscellaneous.name + '.bat'):
                os.remove(md.miscellaneous.name + '.bat')

    return md
    
def _process_love_kernels(field, md):
    """
    This function converts Love kernel solutions from the ISSM (Ice Sheet System Model)
    internal format to a standardized format suitable for further analysis. The conversion
    involves rescaling the kernel components using reference values and applying appropriate
    boundary conditions at the surface layer.

    Parameters
    ----------
    field : numpy.ndarray
        Input Love kernel field data in ISSM format with shape (M, N) where M is the
        flattened spatial dimension and N is the frequency dimension.
    md : object
        Model data structure containing material properties and Love number configuration.
        Must have the following attributes:
        - materials.numlayers : int
            Number of material layers
        - materials.radius : numpy.ndarray
            Radial coordinates of layer boundaries
        - materials.density : numpy.ndarray
            Density values for each layer
        - love.sh_nmax : int
            Maximum spherical harmonic degree
        - love.nfreq : int
            Number of frequencies
        - love.r0 : float
            Reference radius scaling factor
        - love.g0 : float
            Reference gravity scaling factor
        - love.mu0 : float
            Reference shear modulus scaling factor
        - love.forcing_type : int
            Type of forcing (9 or 11 for different boundary conditions)

    Returns
    -------
    numpy.ndarray
        Processed Love kernel field with shape (degmax+1, nfreq, nlayer+1, 6).
        The last dimension contains 6 kernel components:
        - Component 0: Scaled displacement kernel (mm=4)
        - Component 1: Scaled stress kernel (mm=5)
        - Component 2: Scaled displacement kernel (mm=6)
        - Component 3: Scaled stress kernel (mm=1)
        - Component 4: Scaled gravitational potential kernel (mm=2)
        - Component 5: Scaled gravity kernel (mm=3)

    Notes
    -----
    The function handles two main regions:
    1. Interior layers (kk < nlayer): All 6 components are computed from the input field
       with appropriate scaling factors.
    2. Surface boundary layer (kk == nlayer): Special boundary conditions are applied
       depending on the forcing type.
    For forcing_type == 9: Surface stress components are set to zero.
    For forcing_type == 11: Surface stress includes average density correction.
    The average density is computed as a volume-weighted average across all layers.
    """

    nlayer = md.materials.numlayers
    degmax = md.love.sh_nmax
    nfreq = md.love.nfreq
    r0 = md.love.r0
    g0 = md.love.g0
    mu0 = md.love.mu0
    rr = md.materials.radius
    rho = md.materials.density
    forcing_type = md.love.forcing_type

    # Average density term
    rho_avg_partial = np.diff(rr**3, n=1, axis=0)
    rho_avg = np.sum((rho * rho_avg_partial) / rho_avg_partial.sum())

    # Initialize output array: (degree, frequency, layer, components)
    temp_field = np.zeros((degmax + 1, nfreq, nlayer + 1, 6), dtype=float)

    for ii in range(degmax + 1):
        for jj in range(nfreq):
            for kk in range(nlayer + 1):

                if kk < nlayer:
                    
                    # Compute base index in flattened array
                    ll = ii * (nlayer + 1) * 6 + (kk * 6 + 1) + 3

                    temp_field[ii, jj, kk, 0] = field[ll + 0, jj] * r0       # mm = 4
                    temp_field[ii, jj, kk, 1] = field[ll + 1, jj] * mu0      # mm = 5
                    temp_field[ii, jj, kk, 2] = field[ll + 2, jj] * r0       # mm = 6
                    temp_field[ii, jj, kk, 3] = field[ll + 3, jj] * mu0      # mm = 1
                    temp_field[ii, jj, kk, 4] = field[ll + 4, jj] * r0 * g0  # mm = 2
                    temp_field[ii, jj, kk, 5] = field[ll + 5, jj] * g0       # mm = 3

                else:
                    # Surface boundary layer
                    ll = (ii + 1) * (nlayer + 1) * 6 - 2

                    temp_field[ii, jj, kk, 0] = field[ll + 0, jj] * r0
                    temp_field[ii, jj, kk, 2] = field[ll + 1, jj] * r0
                    temp_field[ii, jj, kk, 4] = field[ll + 2, jj] * r0 * g0
                    temp_field[ii, jj, kk, 3] = 0.0  # surface BC

                    if forcing_type == 9:
                        temp_field[ii, jj, kk, 1] = 0.0
                        temp_field[ii, jj, kk, 5] = ((2 * (ii + 1) - 1) / r0
                            - (ii + 1) * field[ll + 2, jj] * g0)

                    elif forcing_type == 11:
                        temp_field[ii, jj, kk, 1] = -(2 * ii + 1) * rho_avg / 3
                        temp_field[ii, jj, kk, 5] = ((2 * (ii + 1) - 1) / r0
                            - (ii + 1) * field[ll + 2, jj] * g0)

        return temp_field

def load_results_from_disk(md, filename):
    """
    Load simulation results from disk into a model object.

    This function loads results from a binary file or processes QMU (Quantification 
    of Margins and Uncertainties) results depending on the model configuration. 
    It handles result parsing, error checking, log file reading, and result 
    organization within the model structure.

    Parameters
    ----------
    md : object
        Model object containing simulation parameters, results structure, and 
        configuration settings. Must have attributes like `qmu`, `results`, 
        `miscellaneous`, and `private`.
    filename : str
        Path to the binary file containing simulation results to be loaded.

    Returns
    -------
    object
        The modified model object with results loaded into the `md.results` 
        attribute and solution type set in `md.private.solution`.

    Raises
    ------
    ValueError
        If no results are found in the specified file.

    Warnings
    --------
    UserWarning
        - If SolutionType field cannot be found in results
        - If errors are detected in the error log during solution

    Notes
    -----
    The function performs different operations based on whether QMU is enabled:
    - **Non-QMU mode**: Parses binary results, loads error/output logs, and 
      organizes results by solution type
    - **QMU mode**: Delegates to QMU-specific postprocessing
    For non-QMU results, the function:
    1. Verifies the result file exists
    2. Initializes the results structure if needed
    3. Parses results and determines solution type
    4. Loads associated log files (.errlog and .outlog)
    5. Extracts single solutions from lists for user convenience
    Log files are read as lists of strings with newlines removed. If no log 
    files exist, empty strings are assigned to the respective log attributes.

    Examples
    --------
    >>> md = load_results_from_disk(md, 'simulation_results.bin')
    >>> print(md.results.StressbalanceSolution.Vel)
    """
    
    if not md.qmu.isdakota:
        # Check that file exists
        if not os.path.exists(filename):
            err_msg = '==========================================================================\n'
            err_msg += '   Binary file {} not found                                              \n'.format(filename)
            err_msg += '                                                                         \n'
            err_msg += '   This typically results from an error encountered during the simulation\n'
            err_msg += '   Please check for error messages above or in the outlog                \n'
            err_msg += '=========================================================================\n'
            print(err_msg)
            return
        
        # Initialize md.results if it is not a structure yet
        if not isinstance(md.results, model.classes.results.default):
            md.results = model.classes.results.default()

        # Load results onto model
        results = parse_results_from_disk(md, filename)
        if not results:
            raise ValueError(f'load_results_from_disk error: no results found in file {filename!r}')
        
        if not hasattr(results[0], 'SolutionType'):
            if hasattr(results[-1], 'SolutionType'):
                results[0].SolutionType = results[-1].SolutionType
            else:
                warnings.warn('load_results_from_disk: Cannot find a SolutionType field in results. Setting to "NoneSolution".')
                results[0].SolutionType = "NoneSolution"
        setattr(md.results, results[0].SolutionType, results)

        # Recover solution_type from results
        md.private.solution = results[0].SolutionType

        # Read log files onto fields
        if os.path.exists(md.miscellaneous.name + '.errlog'):
            with open(md.miscellaneous.name + '.errlog', 'r') as f:
                setattr(getattr(md.results, results[0].SolutionType)[0], 'errlog', [line[:-1] for line in f])
        else:
            setattr(getattr(md.results, results[0].SolutionType)[0], 'errlog', '')

        if os.path.exists(md.miscellaneous.name + '.outlog'):
            with open(md.miscellaneous.name + '.outlog', 'r') as f:
                setattr(getattr(md.results, results[0].SolutionType)[0], 'outlog', [line[:-1] for line in f])
        else:
            setattr(getattr(md.results, results[0].SolutionType)[0], 'outlog', '')

        if getattr(md.results, results[0].SolutionType)[0].errlog:
            warnings.warn('load_results_from_disk: error during solution. Check your errlog and outlog model fields.')

        # If only one solution, extract it from list for user friendliness
        if len(results.steps) == 1 and results[0].SolutionType != 'TransientSolution':
            setattr(md.results, results[0].SolutionType, results[0])

    else:
        # Postprocess QMU results
        md = postprocess_qmu(md, filename)

    return md


def parse_results_from_disk(md, filename):
    """
    Parse and load results from an ISSM binary file.

    This function reads an ISSM binary results file and organizes all fields by
    solution step. It handles sequential reading of field records, groups them
    by time step, and returns a structured results object containing all solution
    data.

    Parameters
    ----------
    md : object
        Model data structure containing simulation parameters and configuration.
        Must have a `constants` attribute with `yts` (years to seconds conversion).
    filename : str
        Path to the binary results file to be parsed.

    Returns
    -------
    model.classes.results.solution
        Structured results object containing all solution steps. Each step
        contains the fields and metadata for that particular solution time.
        Fields are organized as attributes accessible via dot notation.

    Raises
    ------
    IOError
        If the specified results file cannot be opened for reading.
    Exception
        If no results are found in the binary file or if there are issues
        reading the binary data format.

    Notes
    -----
    The function performs the following operations:
    1. Opens and sequentially reads all field records from the binary file
    2. Extracts unique solution steps and organizes fields by step
    3. Applies unit conversions and data formatting (e.g., flattening column vectors)
    4. Returns a results object with fields accessible by step index

    The binary file format is expected to contain field records with names,
    time stamps, step numbers, and data arrays in ISSM-specific format.
    """

    try:
        fid = open(filename, 'rb')
    except IOError:
        raise IOError(f"parse_results_from_disk: could not open {filename}")

    # --- read all fields sequentially ---
    all_fields = []
    while True:
        field_record = _read_data(fid, md)
        if field_record is None:
            break
        all_fields.append(field_record)
    fid.close()

    if not all_fields:
        raise Exception(f'parse_results_from_disk: no results found in binary file {filename}')

    # --- determine all unique steps ---
    all_steps = np.array([r['step'] for r in all_fields if r['step'] != -9999])
    all_steps = np.sort(np.unique(all_steps))

    # --- construct solution array ---
    results = model.classes.results.solution()
    
    # Assign fields to appropriate solution step
    for i in range(len(all_fields)):
        result = all_fields[i]
        index = 0
        if result['step'] != -9999:
            index = np.where(result['step'] == all_steps)[0][0]
            setattr(results[index], 'step', result['step'])
        if result['time'] != -9999:
            setattr(results[index], 'time', result['time'])

        field = result['field']

        # Flatten 2D arrays for improved python indexing
        if isinstance(field, np.ndarray) and field.ndim == 2 and field.shape[1] == 1:
            # Convert column vectors to 1D arrays
            setattr(results[index], result['fieldname'], field.ravel())
        else:
            setattr(results[index], result['fieldname'], field)

    return results


def _read_data(fid, md):
    """
    Read a single field record from an ISSM binary results file.

    This function reads one complete field record from an ISSM (Ice Sheet System Model)
    binary results file in little-endian format. Each record contains metadata
    (field name, time, step, data type) followed by the actual field data.

    Parameters
    ----------
    fid : file object
        File handle opened in binary read mode pointing to an ISSM results file.
    md : object
        Model data structure containing constants and configuration parameters.
        Must have a `constants.yts` attribute for time unit conversion.

    Returns
    -------
    collections.OrderedDict or None
        Dictionary containing the parsed field record with keys:
        - 'fieldname' : str
            Name of the field being read
        - 'time' : float
            Time value associated with the field (converted from seconds to years)
        - 'step' : int
            Solution step number
        - 'field' : numpy.ndarray or str
            The actual field data, with appropriate unit conversions applied
        Returns None if end-of-file is reached or if an error occurs during reading.

    Notes
    -----
    The function handles multiple data types:
    - Type 1: Double precision vector
    - Type 2: String data
    - Type 3: Double precision matrix
    - Type 4: Integer matrix
    
    Unit conversions are automatically applied based on field names using the
    `_convert_units` helper function. Time values are normalized from seconds
    to years using the years-to-seconds conversion factor from the model.

    The binary format expects little-endian byte ordering and reads data in
    the following sequence:
    1. Field name length (4 bytes, int32)
    2. Field name (variable length string)
    3. Time value (8 bytes, float64)
    4. Step number (4 bytes, int32)
    5. Data type (4 bytes, int32)
    6. Data dimensions and values (format depends on data type)
    """

    try:
        # --- read field name length ---
        length_bytes = np.fromfile(fid, dtype='<i4', count=1)
        if length_bytes.size == 0:
            return None  # EOF
        length = int(length_bytes[0])

        # --- read field name ---
        fieldname_bytes = np.fromfile(fid, dtype='S1', count=length)
        fieldname = b''.join(fieldname_bytes).rstrip(b'\x00').decode()

        # --- read metadata ---
        time = np.fromfile(fid, dtype='<f8', count=1)[0]
        step = np.fromfile(fid, dtype='<i4', count=1)[0]
        datatype = np.fromfile(fid, dtype='<i4', count=1)[0]
        M = np.fromfile(fid, dtype='<i4', count=1)[0]

        # --- read data ---
        if datatype == 1:  # double vector
            field = np.fromfile(fid, dtype='<f8', count=M)
        elif datatype == 2:  # string
            field_bytes = np.fromfile(fid, dtype='S1', count=M)
            field = b''.join(field_bytes).rstrip(b'\x00').decode()
        elif datatype == 3:  # double matrix
            N = int(np.fromfile(fid, dtype='<i4', count=1)[0])
            field = np.fromfile(fid, dtype='<f8', count=M*N).reshape(M, N)
        elif datatype == 4:  # int matrix
            N = int(np.fromfile(fid, dtype='<i4', count=1)[0])
            field = np.fromfile(fid, dtype='<i4', count=M*N).reshape(M, N)
        else:
            raise TypeError(f"cannot read data of datatype {datatype}")

        # --- unit conversions ---
        field = _convert_units(fieldname, field, md)

        # --- normalize time ---
        if time != -9999:
            time /= md.constants.yts

        return collections.OrderedDict(fieldname=fieldname, time=time, step=step, field=field)

    except Exception:
        # any read error or EOF mid-field
        return None
    
def _convert_units(field_name, field, md):
    """
    Apply unit conversions based on field name.

    This function converts ISSM simulation output fields from their internal
    units to more convenient units for analysis. It handles velocity fields,
    mass balance terms, and specialized fields like Love kernels.

    Parameters
    ----------
    field_name : str
        Name of the field being converted. Used to determine which conversion
        factors to apply.
    field : numpy.ndarray or scalar
        The field data to be converted. Can be any numeric type that supports
        arithmetic operations.
    md : object
        Model data structure containing conversion constants. Must have a
        `constants.yts` attribute (years to seconds conversion factor).

    Returns
    -------
    numpy.ndarray or scalar
        Field data with appropriate unit conversions applied:
        - Velocity fields: converted from m/s to m/yr
        - Mass balance fields: converted from kg/s to m/yr equivalent
        - Total/cumulative fields: converted to Gt/yr
        - Love kernels: processed using specialized kernel transformation
        - Other fields: returned unchanged

    Notes
    -----
    Unit conversions applied:
    - **Time-based fields** (velocities, rates): multiplied by yts to convert
      from per-second to per-year units
    - **Total/cumulative fields**: multiplied by yts/1e12 to convert from kg/s
      to Gt/yr (gigatons per year)
    - **Love kernels**: processed using `_process_love_kernels` for specialized
      geophysical unit handling

    The following field categories are recognized:
    - Velocity components (Vx, Vy, Vz, Vel, etc.)
    - Surface mass balance terms (SmbMassBalance, SmbPrecipitation, etc.)
    - Calving rates (CalvingCalvingrate, Calvingratex, etc.)
    - Total mass fluxes (TotalFloatingBmb, GroundinglineMassFlux, etc.)

    Examples
    --------
    >>> # Convert velocity from m/s to m/yr
    >>> vel_converted = _convert_units('Vx', velocity_field, md)
    >>> # Convert mass balance from kg/s to Gt/yr  
    >>> mass_converted = _convert_units('TotalSmb', mass_field, md)
    """
    yts = md.constants.yts
    gigaton_scale = yts / 1e12

    # Time-based fields
    scale_fields = {
        "BalancethicknessThickeningRate", "HydrologyWaterVx", "HydrologyWaterVy",
        "Vx", "Vy", "Vz", "Vel", "VxShear", "VyShear",
        "VxBase", "VyBase", "VxSurface", "VySurface",
        "VxAverage", "VyAverage", "VxDebris", "VyDebris",
        "SmbMassBalance", "SmbPrecipitation", "SmbRain", "SmbRunoff", "SmbRunoffSubstep",
        "SmbEvaporation", "SmbRefreeze", "SmbEC", "SmbAccumulation", "SmbMelt",
        "SmbMAdd", "SmbWAdd", "CalvingCalvingrate", "Calvingratex", "Calvingratey",
        "CalvingMeltingrate",
    }

    # Cumulative / Total Fields (Gt yr^{-1})
    total_fields = {
        "TotalFloatingBmb", "TotalFloatingBmbScaled",
        "TotalGroundedBmb", "TotalGroundedBmbScaled",
        "TotalSmb", "TotalSmbScaled",
        "TotalSmbMelt", "TotalSmbRefreeze",
        "GroundinglineMassFlux", "IcefrontMassFlux", "IcefrontMassFluxLevelset",
    }

    if field_name in scale_fields:
        return field * yts
    elif field_name in total_fields:
        return field * gigaton_scale
    elif field_name in ['LoveKernelsReal', 'LoveKernelsImag']:
        field = _process_love_kernels(field_name, field, md)
    else:
        return field
