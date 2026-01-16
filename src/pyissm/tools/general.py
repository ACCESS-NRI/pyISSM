"""
Utility functions for ISSM

This module contains various utility functions that are used throughout the ISSM codebase.
"""

import numpy as np
import struct
import math
from pyissm import model

## ------------------------------------------------------------------------------------
## UNIT CONVERSIONS
## ------------------------------------------------------------------------------------
def convert_units(input_units,
                  output_units,
                  data,
                  yts = 365 * 24 * 60 * 60,
                  rho_ice = 917):
    """
    Convert numerical data between supported geophysical units.

    This function supports a variety of conversions commonly required when
    working with ISSM input/output entities. Required constants are
    consistent with those used in the Ice sheet and Sea Level System Model (ISSM).

    Parameters
    ----------
    input_units : str
        Units of the input data. Must be one of:
        'm', 'km', 'ms-1', 'myr-1', 'm2', 'km2', 'Gt', 'km3', 'Gtyr-1', 'kgs-1'.
    output_units : str
        Desired units for the output data. Must be one of:
        'm', 'km', 'ms-1', 'myr-1', 'm2', 'km2', 'Gt', 'km3', 'Gtyr-1', 'kgs-1'.
    data : float, list, or ndarray
        Numerical value(s) to convert. Data are converted to array for computation.
    yts : float or int, optional
        Seconds in a year (default = 365 * 24 * 60 * 60)
    rho_ice : float or int, optional
        Ice density in kg/m3 (default: 917)

    Returns
    -------
    converted_data : float or ndarray
        The data converted from `input_units` to `output_units`.

    Raises
    ------
    ValueError
        If `input_units` or `output_units` are not among the accepted units,
        or if the requested unit conversion is not supported.

    Notes
    -----
    The following unit conversions are supported:

    - Length: m <-> km
    - Area: m² <-> km²
    - Speed: m/s <-> m/yr
    - Volume: m³ <-> km³
    - Mass: Gt <-> km³ (using ice density 917 kg/m³)
    - Rate: Gt/yr <-> kg/s
    - Rate: kg/m²/s-1 <-> myr-1
    - Rate: kg/m²/s-1 <-> myr-1ie
    """

    ## Define list of units supported by function
    accepted_units = ['m', 'km',
                      'ms-1', 'myr-1',
                      'm2', 'km2',
                      'Gt', 'km3',
                      'm3', 'kg',
                      'Gtyr-1', 'kgs-1',
                      'myr-1ie', 'myr-1',
                      'kgm-2s-1']

    ## Argument error checks
    if input_units not in accepted_units:
        raise ValueError(f"convert_units: Invalid input_units '{input_units}'. Must be one of {accepted_units}")

    if output_units not in accepted_units:
        raise ValueError(f"convert_units: Invalid output_units '{output_units}'. Must be one of {accepted_units}")

    ## Convert to array for element-wise math
    data = np.asarray(data, dtype = np.float64)

    ## Length conversions -----------------------
    if input_units == 'm' and output_units == 'km':
        converted_data = data / 1e3

    elif input_units == 'km' and output_units == 'm':
        converted_data = data * 1e3

    ## Area conversions -----------------------
    elif input_units == 'm2' and output_units == 'km2':
        converted_data = data / 1e6

    elif input_units == 'km2' and output_units == 'm2':
        converted_data = data * 1e6

    ## Speed conversions -----------------------
    elif input_units == 'ms-1' and output_units == 'myr-1':
        converted_data = data * yts

    elif input_units == 'myr-1' and output_units == 'ms-1':
        converted_data = data / yts

    ## Volume conversions -----------------------
    elif input_units == 'm3' and output_units == 'km3':
        converted_data = data / 1e9

    elif input_units == 'km3' and output_units == 'm3':
        converted_data = data * 1e9

    ## Mass conversions -----------------------
    elif input_units == 'Gt' and output_units == 'km3':
        kg = data * 1e12
        m3 = kg / rho_ice
        km3 = m3 / 1e9
        converted_data = km3

    elif input_units == 'km3' and output_units == 'Gt':
        m3 = data * 1e9
        kg = m3 * rho_ice
        gt = kg / 1e12
        converted_data = gt

    elif input_units == 'm3' and output_units == 'kg':
        converted_data = data * rho_ice

    elif input_units == 'kg' and output_units == 'm3':
        converted_data = data / rho_ice

    ## Rate conversions -----------------------
    elif input_units == 'Gtyr-1' and output_units == 'kgs-1':
        converted_data = (data * 1e12) / yts

    elif input_units == 'kgs-1' and output_units == 'Gtyr-1':
        converted_data = (data * yts) / 1e12

    elif input_units == 'myr-1ie' and output_units == 'kgm-2s-1':
        converted_data = data / yts

    elif input_units == 'kgm-2s-1' and output_units == 'myr-1ie':
        converted_data = data * yts

    elif input_units == 'myr-1' and output_units == 'kgm-2s-1':
        converted_data = (data / yts) * rho_ice

    elif input_units == 'kgm-2s-1' and output_units == 'myr-1':
        converted_data = (data / rho_ice) * yts

    else:
        raise ValueError("""convert_units: Requested unit conversion currently not supported. Available conversions include:
- m <-> km
- m2 <-> km2
- ms-1 <-> my-1
- m3 <-> km3
- Gt <-> km3
- m3 <-> kg
- Gtyr-1 <-> kgs-1
- myr-1 <-> kgm-2s-1
- myr-1ie <-> kgm-2s-1""")

    return converted_data

## ------------------------------------------------------------------------------------
## ISSM MODEL UTILITIES
## ------------------------------------------------------------------------------------
def has_nested_attr(obj, *attrs):
    """
    Check whether an object has a chain of nested attributes.

    This function tests whether an object has a sequence of attributes,
    each accessible from the previous one, such as `obj.a.b.c`.

    Parameters
    ----------
    obj : object
        The base object from which to start the attribute lookup.
    *attrs : str
        A sequence of attribute names representing the path to check.
        Each string in `attrs` should be a valid attribute of the previous one.

    Returns
    -------
    has_attr : bool
        `True` if all nested attributes exist, `False` otherwise.
    """
    for attr in attrs:
        if not hasattr(obj, attr):
            return False
        obj = getattr(obj, attr)
    return True
  
def extract_field_layer(md,
                        field,
                        layer):
    """
    Extract a 2D horizontal layer from a 3D field defined on the model's vertices.

    This function isolates a specific horizontal layer from a field defined over all
    3D vertices of a model. Only vertex-based fields are supported; element-based
    fields will raise an error.

    Parameters
    ----------
    md : ISSM Model object
        ISSM Model object containing mesh. Must be compatible with process_mesh().
    field : np.ndarray
        A 1D array of shape 'md.mesh.numberofvertices' representing a 3D field defined on vertices.
    layer : int
        The vertical layer index to extract (1-based). Must be a single integer.
        Layer indexing starts from 1 (base) to 'md.mesh.numberoflayers' (surface).

    Returns
    -------
    select_layer : np.ndarray
        A 1D array of shape 'md.mesh.numberofvertices2d' or 'md.mesh.numberofelements2d' containing the field values at the specified layer.
    select_indices : np.ndarray
        A 1D array of shape 'md.mesh.numberofvertices2d' or 'md.mesh.numberofelements2d' containing the indices of vertices/elements associated with the specified layer.

    Notes
    -----
    - This function assumes that the 3D vertex/element ordering in 'field' is layer-major:
      i.e., all vertices/elements in layer 1 come first, then layer 2, and so on.
    - Depth-averaging functionality (e.g., specifying a range of layers) is not yet implemented.

    Example
    -----
    field_2d, select_indices = extract_field_layer(md, md.results.TransientSolution.Temperature[0], layer = 1)
    field_2d, select_indices = extract_field_layer(md, md.materials.rheology_n, layer = 6)
    """

    ## Error Checks
    # Process mesh
    mesh, mesh_x, mesh_y, mesh_elements, is3d = model.mesh.process_mesh(md)

    # If it's not a 3D model, throw and error
    if not is3d:
        raise TypeError('extract_field_layer: provided model is not 3D')

    # TODO: Implement "layer = [0, 8]" to DepthAverage over a given layer range.
    # Need to 'weight' the average based on the layer thickness
    if isinstance(layer, list):
        raise ValueError('extract_field_layer: A single numeric layer must be defined. Depth averaging is not yet supported.')

    # Check dimensions of model mesh. The number of vertices must be equal to md.mesh.numberofvertices2d * md.mesh.numberoflayers
    if not np.equal(md.mesh.numberofvertices / md.mesh.numberoflayers, md.mesh.numberofvertices2d):
        raise ValueError('extract_field_layer: model mesh is not correctly dimensioned')

    ## Extract data defined on elements
    if field.shape[0] == md.mesh.numberofelements:

        # The number of element layers is (md.mesh.numberoflayers - 1). Check that the defined layer can be extracted.
        if not (1 <= layer <= md.mesh.numberoflayers - 1):
            raise ValueError(f'Layer {layer} is out of bounds for element-based fields. Select a layer between 1 and {md.mesh.numberoflayers -1}')

        # Identify element positions for requested layer
        start_pos = md.mesh.numberofelements2d * (layer - 1)
        end_pos = md.mesh.numberofelements2d * layer
        select_indices = np.arange(start_pos, end_pos)

        # Extract requested layer from field
        select_layer = field[select_indices]

    ## Extract data defined on vertices
    if field.shape[0] == md.mesh.numberofvertices:
        # Identify vertex positions for requested layer
        start_pos = md.mesh.numberofvertices2d * (layer - 1)
        end_pos = md.mesh.numberofvertices2d * layer
        select_indices = np.arange(start_pos, end_pos)

        ## Extract requested layer from field
        select_layer = field[select_indices]

    return select_layer, select_indices

## ------------------------------------------------------------------------------------
## SOLID EARTH FUNCTION
## ------------------------------------------------------------------------------------
def planetradius(planet):
    """
        Return the mean radius of a specified planetary body.

        Parameters
        ----------
        planet : str
            Name of the planet. Supported values are:
            - `'earth'` : Earth's mean radius in meters.
            - `'europa'` : Europa's mean radius in meters.

        Returns
        -------
        radius : float
            Planetary radius in meters.

        Raises
        ------
        TypeError
            If the input `planet` is not one of the supported values.
    """

    if planet == 'earth':
        radius = 6.371012e6
    elif planet == 'europa':
        radius = 1.5008e6
    else:
        raise TypeError(f'planetradius: Planet type {planet} is not supported.')

    return radius

def _wgs84_ellipsoid_constants():
    """
    Return constants for the WGS84 ellipsoid used in coordinate conversions.

    Returns
    -------
    re : float
        Equatorial radius of the WGS84 ellipsoid in meters.
    f : float
        Flattening of the WGS84 ellipsoid.
    ex2 : float
        Eccentricity squared of the WGS84 ellipsoid.
    ex : float
        Eccentricity of the WGS84 ellipsoid.
    """
    
    # Equitorial radius of the earth in meters
    re = 6378137

    # Flattening for WGS84 ellipsoid
    f  = 1./298.257223563
    
    # Eccentricity squared
    ex2 = 2*f - f**2

    # Eccentricity
    ex = np.sqrt(ex2)

    return re, f, ex2, ex

def xy_to_ll(x,
             y,
             sign,
             central_meridian = None,
             standard_parallel = None):
    
    """
    Convert Cartesian coordinates (x, y) to geographic coordinates (latitude, longitude)
    using the polar stereographic projection.

    Parameters
    ----------
    x : array_like
        The x-coordinates in meters.
    y : array_like
        The y-coordinates in meters.
    sign : int
        The hemisphere indicator, where 1 indicates the northern hemisphere and -1 indicates the southern hemisphere.
    central_meridian : float, optional
        The central meridian in degrees. If not specified, defaults are used based on the hemisphere.
    standard_parallel : float, optional
        The standard parallel in degrees. If not specified, defaults are used based on the hemisphere.

    Returns
    -------
    lat : ndarray
        The latitude in degrees.
    lon : ndarray
        The longitude in degrees, adjusted by the central meridian.

    Raises
    ------
    ValueError
        If `sign` is not 1 or -1.
        If only one of `central_meridian` or `standard_parallel` is specified.

    Notes
    -----
    The function uses the WGS84 ellipsoid parameters for the conversion.
    Special handling is included for the case when `standard_parallel` is 90 degrees.
    """

    # Error checks
    if sign not in [1, -1]:
        raise ValueError('pyissm.tools.general.xy_to_ll: sign should be either 1 (north) or -1 (south)')
    
    # Set defaults depending on hemisphere
    if central_meridian is None and standard_parallel is None:
        if sign == 1:
            central_meridian = 45.
            standard_parallel = 70.
            print('pyissm.tools.general.xy_to_ll: creating coordinates in north polar stereographic (Std Latitude: 70degN Meridian: 45deg)')

        elif sign == -1:
            central_meridian = 0.
            standard_parallel = 71.
            print('pyissm.tools.general.xy_to_ll: creating coordinates in south polar stereographic (Std Latitude: 71degS Meridian: 0deg)')

    elif (central_meridian is None) != (standard_parallel is None):
        raise ValueError("Specify both central_meridian and standard_parallel, or neither.")
    
    # Ensure x and y are numpy arrays
    x = np.asarray(x, dtype = float)
    y = np.asarray(y, dtype = float)

    # Define constants
    re, _, ex2, ex = _wgs84_ellipsoid_constants()

    # Convert
    sl = np.deg2rad(standard_parallel)
    rho = np.hypot(x, y)

    cm = np.cos(sl) / np.sqrt(1 - ex2 * np.sin(sl)**2)
    T = np.tan(np.pi/4 - sl/2) / ((1 - ex*np.sin(sl)) / (1 + ex*np.sin(sl)))**(ex/2)

    # Special case: standard_parallel = 90deg
    if np.isclose(standard_parallel, 90.0):
        T = rho * np.sqrt((1 + ex)**(1 + ex) * (1 - ex)**(1 - ex)) / (2 * re)
    else:
        T = rho * T / (re * cm)

    chi = np.pi/2 - 2 * np.arctan(T)

    lat = (
        chi
        + (ex2/2 + 5*ex2**2/24 + ex2**3/12) * np.sin(2*chi)
        + (7*ex2**2/48 + 29*ex2**3/240) * np.sin(4*chi)
        + (7*ex2**3/120) * np.sin(6*chi)
    )

    lat *= sign
    lon = sign * np.arctan2(sign*x, -sign*y)

    # Handle near-origin
    near_origin = rho <= 0.1
    if np.any(near_origin):
        lat = lat.copy()
        lon = lon.copy()
        lat[near_origin] = sign * (np.pi/2)
        lon[near_origin] = 0.0

    lat = np.rad2deg(lat)
    lon = np.rad2deg(lon) - central_meridian

    return lat, lon


def ll_to_xy(lat,
             lon,
             sign,
             central_meridian = None,
             standard_parallel = None):
    """
    Convert geographic coordinates (latitude, longitude) to Cartesian coordinates (x, y)
    using the polar stereographic projection.

    Parameters
    ----------
    lat : array_like
        The latitude in degrees.
    lon : array_like
        The longitude in degrees.
    sign : int
        The hemisphere indicator, where 1 indicates the northern hemisphere and -1 indicates the southern hemisphere.
    central_meridian : float, optional
        The central meridian in degrees. If not specified, defaults are used based on the hemisphere.
    standard_parallel : float, optional
        The standard parallel in degrees. If not specified, defaults are used based on the hemisphere.

    Returns
    -------
    x : ndarray
        The x-coordinates in meters.
    y : ndarray
        The y-coordinates in meters.

    Raises
    ------
    ValueError
        If `sign` is not 1 or -1.
        If only one of `central_meridian` or `standard_parallel` is specified.

    Notes
    -----
    The function uses the WGS84 ellipsoid parameters for the conversion.
    Special handling is included for the case when `standard_parallel` is 90 degrees.
    """

    # Error checks
    if sign not in [1, -1]:
        raise ValueError('ll_to_xy: sign should be either 1 (north) or -1 (south)')

    # Set defaults depending on hemisphere
    if central_meridian is None and standard_parallel is None:
        if sign == 1:
            central_meridian = 45.
            standard_parallel = 70.
            print('ll_to_xy: using north polar stereographic (Std Lat: 70N, Meridian: 45E)')
        else:
            central_meridian = 0.
            standard_parallel = 71.
            print('ll_to_xy: using south polar stereographic (Std Lat: 71S, Meridian: 0E)')
    elif (central_meridian is None) != (standard_parallel is None):
        raise ValueError("Specify both central_meridian and standard_parallel, or neither.")

    # Ensure lat/lon are numpy arrays
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)

    # Define constants
    re, _, ex2, ex = _wgs84_ellipsoid_constants()

    # Convert degrees to radians
    lat_rad = np.deg2rad(np.abs(lat))
    lon_rad = np.deg2rad(lon + central_meridian)

    # Compute T
    T = np.tan(np.pi/4 - lat_rad/2) / ((1 - ex*np.sin(lat_rad)) / (1 + ex*np.sin(lat_rad)))**(ex/2)

    # Standard parallel calculations
    sl_rad = np.deg2rad(standard_parallel)
    if np.isclose(standard_parallel, 90.0):
        rho = 2 * re * T / np.sqrt((1 + ex)**(1 + ex) * (1 - ex)**(1 - ex))
    else:
        mc = np.cos(sl_rad) / np.sqrt(1 - ex2 * np.sin(sl_rad)**2)
        Tc = np.tan(np.pi/4 - sl_rad/2) / ((1 - ex*np.sin(sl_rad)) / (1 + ex*np.sin(sl_rad)))**(ex/2)
        rho = re * mc * T / Tc

    # Convert to X, Y
    x = rho * sign * np.sin(sign * lon_rad)
    y = -rho * sign * np.cos(sign * lon_rad)

    # Handle near-pole points
    near_pole = np.abs(lat_rad - np.pi/2) < 1e-10
    if np.any(near_pole):
        x = x.copy()
        y = y.copy()
        x[near_pole] = 0.0
        y[near_pole] = 0.0

    return x, y

def compare_binfiles(file1, file2, verbose=0, outfile=None):
    """
    Compare two ISSM-style binary files field by field.

    Parameters
    ----------
    file1 : str
        Path to the first binary file.
    file2 : str
        Path to the second binary file.
    verbose : int, optional
        Verbosity level (2 = show differing array indices).
    outfile : str or None, optional
        Path to write the output table. If None, prints to terminal.

    Returns
    -------
    rows : list of tuples
        List of (field, status, file1_summary, file2_summary)
    """

    # ------------------------------
    # Internal helpers
    # ------------------------------
    def _code_to_format(code):
        mapping = {
            1:'Boolean', 2:'Integer', 3:'Double', 4:'String', 5:'BooleanMat',
            6:'IntMat', 7:'DoubleMat', 8:'MatArray', 9:'StringArray'
        }
        return mapping.get(code, None)

    def _read_bin_to_dict(infile):
        data_dict = {}
        with open(infile, 'rb') as f:
            while True:
                try:
                    recordnamesize = struct.unpack('i', f.read(4))[0]
                except struct.error:
                    break  # EOF

                recordname = struct.unpack(f'{recordnamesize}s', f.read(recordnamesize))[0].decode('ASCII')
                reclen = struct.unpack('q', f.read(8))[0]
                code = struct.unpack('i', f.read(4))[0]
                fmt = _code_to_format(code)

                if fmt in ['Boolean', 'Integer']:
                    val = struct.unpack('i', f.read(reclen - 4))[0]
                elif fmt == 'Double':
                    val = struct.unpack('d', f.read(reclen - 4))[0]
                elif fmt == 'String':
                    strlen = struct.unpack('i', f.read(4))[0]
                    val = struct.unpack(f'{strlen}s', f.read(strlen))[0].decode('ASCII')
                elif fmt in ['BooleanMat', 'IntMat', 'DoubleMat']:
                    _ = struct.unpack('i', f.read(4))[0]  # mattype
                    s0 = struct.unpack('i', f.read(4))[0]
                    s1 = struct.unpack('i', f.read(4))[0]
                    mat = np.zeros((s0, s1))
                    for i in range(s0):
                        for j in range(s1):
                            mat[i, j] = struct.unpack('d', f.read(8))[0]
                    val = mat
                elif fmt in ['MatArray', 'StringArray']:
                    f.seek(reclen - 4, 1)
                    val = None
                else:
                    raise TypeError(f'Unsupported type code {code} ({recordname})')

                data_dict[recordname] = val

        return data_dict

    def _summarize_value(val):
        if val is None:
            return "N/A"

        if isinstance(val, np.ndarray):
            if val.size == 0:
                return f'empty array shape={val.shape}, dtype={val.dtype}'
            if np.issubdtype(val.dtype, np.floating):
                n_nan = np.isnan(val).sum()
                if n_nan == val.size:
                    return f'all-NaN array shape={val.shape}, dtype={val.dtype}, NaNs={n_nan}'
                return (f'array shape={val.shape}, dtype={val.dtype}, '
                        f'min={np.nanmin(val):.3g}, max={np.nanmax(val):.3g}, NaNs={n_nan}')
            else:
                n_nan = 0
                return f'array shape={val.shape}, dtype={val.dtype}, NaNs={n_nan}'

        if isinstance(val, float):
            if math.isnan(val):
                return "NaN"
            return f'{val}'

        return str(val)

    def _compare_values(val1, val2):
        # Arrays
        if isinstance(val1, np.ndarray) and isinstance(val2, np.ndarray):
            if val1.shape != val2.shape:
                return "DIFFERENT", f"shape {val1.shape} vs {val2.shape}"

            both_nan = np.isnan(val1) & np.isnan(val2)
            diff_mask = ~both_nan & (val1 != val2)
            n_diff = np.count_nonzero(diff_mask)

            if n_diff == 0:
                return "SAME", ""
            else:
                detail = f"{n_diff} elements differ"
                if verbose > 1:
                    detail += f" at indices {np.argwhere(diff_mask).tolist()}"
                return "DIFFERENT", detail

        # Scalars
        if isinstance(val1, float) and isinstance(val2, float):
            if math.isnan(val1) and math.isnan(val2):
                return "SAME", ""

        if val1 == val2:
            return "SAME", ""

        return "DIFFERENT", f"{val1} != {val2}"

    # ------------------------------
    # Read files
    # ------------------------------
    dict1 = _read_bin_to_dict(file1)
    dict2 = _read_bin_to_dict(file2)

    all_keys = sorted(set(dict1) | set(dict2))

    # ------------------------------
    # Prepare rows
    # ------------------------------
    rows = []
    for key in all_keys:
        v1 = dict1.get(key, None)
        v2 = dict2.get(key, None)

        if key not in dict1:
            status = "MISSING"
            s1 = "<missing>"
            s2 = _summarize_value(v2)
        elif key not in dict2:
            status = "MISSING"
            s1 = _summarize_value(v1)
            s2 = "<missing>"
        else:
            status, _ = _compare_values(v1, v2)
            s1 = _summarize_value(v1)
            s2 = _summarize_value(v2)

        rows.append((key, status, s1, s2))

    # ------------------------------
    # Column widths for printing
    # ------------------------------
    col_widths = [int(max(len(str(r[i])) for r in rows + [("Field","Status","File1 summary","File2 summary")]))
                  for i in range(4)]

    # ------------------------------
    # Prepare table lines
    # ------------------------------
    output_lines = []
    output_lines.append(f'Comparing:\n  {file1}\n  {file2}\n')
    header = "{:<{w0}}  {:<{w1}}  {:<{w2}}  {:<{w3}}".format(
        "Field", "Status", "File1 summary", "File2 summary",
        w0=col_widths[0], w1=col_widths[1], w2=col_widths[2], w3=col_widths[3]
    )
    output_lines.append(header)
    output_lines.append('-' * (sum(col_widths) + 6))

    for r in rows:
        output_lines.append("{:<{w0}}  {:<{w1}}  {:<{w2}}  {:<{w3}}".format(
            r[0], r[1], r[2], r[3],
            w0=col_widths[0], w1=col_widths[1], w2=col_widths[2], w3=col_widths[3]
        ))

    # ------------------------------
    # Output to terminal or file
    # ------------------------------
    if outfile:
        with open(outfile, 'w') as f:
            for line in output_lines:
                f.write(line + '\n')
    else:
        for line in output_lines:
            print(line)