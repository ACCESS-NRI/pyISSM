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

def compare_bin_files(
    file1,
    file2,
    out_file = None,
    *,
    compare_data = False,
    compare_shape = False,
    compare_format = False,
    compare_mattype = False,
):
    """
    Compare two ISSM binary files.

    This function compares two binary files in ISSM format and reports differences
    based on the selected comparison mode. Exactly one of the comparison modes must
    be enabled.

    Parameters
    ----------
    file1 : str
        Path to the first binary file to compare.
    file2 : str
        Path to the second binary file to compare.
    out_file : str, optional
        Path to write the comparison report. If not provided, output is printed
        to the console. Default is None.
    compare_data : bool, optional
        If True, compare the actual data values in the files. Default is False.
    compare_shape : bool, optional
        If True, compare the shape of array fields. Default is False.
    compare_format : bool, optional
        If True, compare the data type format of fields. Default is False.
    compare_mattype : bool, optional
        If True, compare the matrix type codes of fields. Default is False.

    Returns
    -------
    None
        Output is either written to a file (if `out_file` is specified) or
        printed to the console.

    Raises
    ------
    ValueError
        If exactly one of the compare_* options is not True.
    TypeError
        If an unsupported type code is encountered in the binary file.
    struct.error
        If there is an error reading the binary file format.

    Notes
    -----
    Exactly ONE of the `compare_data`, `compare_shape`, `compare_format`, or
    `compare_mattype` options must be True. An error will be raised if zero or
    more than one option is enabled.

    The comparison report displays fields that are missing in either file,
    along with a status indicating whether values are the same or different
    based on the selected comparison mode.
    """

    # Parse arguments -- enforce one option
    modes = {
        "data": compare_data,
        "shape": compare_shape,
        "format": compare_format,
        "mattype": compare_mattype,
    }
    active = [k for k, v in modes.items() if v]

    if len(active) != 1:
        raise ValueError(
            "Exactly one comparison mode must be True: "
            "compare_data, compare_shape, compare_format, compare_mattype"
        )

    mode = active[0]

    # Define internal helper functions
    ## Map binary file codes to formats
    def _code_to_format(code):
        return {
            1: 'Boolean',
            2: 'Integer',
            3: 'Double',
            4: 'String',
            5: 'BooleanMat',
            6: 'IntMat',
            7: 'DoubleMat',
            8: 'MatArray',
            9: 'StringArray',
        }.get(code)

    ## Unpack binary file to a dictionary
    def _read_bin_to_dict(infile):
        data = {}

        with open(infile, "rb") as f:
            while True:
                try:
                    namesize = struct.unpack("i", f.read(4))[0]
                except struct.error:
                    break

                name = struct.unpack(
                    f"{namesize}s", f.read(namesize)
                )[0].decode("ASCII")

                reclen = struct.unpack("q", f.read(8))[0]
                code = struct.unpack("i", f.read(4))[0]
                fmt = _code_to_format(code)

                mattype = None
                shape = None

                if fmt in ("Boolean", "Integer"):
                    val = struct.unpack("i", f.read(reclen - 4))[0]

                elif fmt == "Double":
                    val = struct.unpack("d", f.read(reclen - 4))[0]

                elif fmt == "String":
                    strlen = struct.unpack("i", f.read(4))[0]
                    val = struct.unpack(
                        f"{strlen}s", f.read(strlen)
                    )[0].decode("ASCII")

                elif fmt in ("BooleanMat", "IntMat", "DoubleMat"):
                    mattype = struct.unpack("i", f.read(4))[0]
                    n0 = struct.unpack("i", f.read(4))[0]
                    n1 = struct.unpack("i", f.read(4))[0]
                    shape = (n0, n1)

                    mat = np.zeros(shape)
                    for i in range(n0):
                        for j in range(n1):
                            mat[i, j] = struct.unpack("d", f.read(8))[0]
                    val = mat

                elif fmt in ("MatArray", "StringArray"):
                    f.seek(reclen - 4, 1)
                    val = None

                else:
                    raise TypeError(f"Unsupported type code {code} ({name})")

                data[name] = dict(
                    value=val,
                    code=code,
                    format=fmt,
                    mattype=mattype,
                    shape=shape,
                )

        return data

    ## Summarise a given field
    def _summarize(record):
        if record is None:
            return "<missing>"

        if mode == "format":
            return f"{record['format']} (code {record['code']})"

        if mode == "mattype":
            return str(record["mattype"])

        if mode == "shape":
            return str(record["shape"])

        # mode == "data"
        val = record["value"]
        if val is None:
            return "N/A"

        if isinstance(val, np.ndarray):
            if val.size == 0:
                return f"empty {val.shape}"
            elif np.all(np.isnan(val)):
                return f"{val.shape}, all NaN"
            return (
                f"{val.shape}, "
                f"min={np.nanmin(val):.3g}, "
                f"max={np.nanmax(val):.3g}, "
                f"NaNs={np.isnan(val).sum()}"
            )

        if isinstance(val, float):
            return "NaN" if math.isnan(val) else str(val)

        return str(val)

    ## Compare two fields
    def _compare(r1, r2):
        if mode == "format":
            return (
                "SAME" if
                (r1["code"], r1["format"]) == (r2["code"], r2["format"])
                else "DIFFERENT"
            )

        if mode == "mattype":
            return "SAME" if r1["mattype"] == r2["mattype"] else "DIFFERENT"

        if mode == "shape":
            return "SAME" if r1["shape"] == r2["shape"] else "DIFFERENT"

        # mode == "data"
        v1, v2 = r1["value"], r2["value"]

        if isinstance(v1, np.ndarray) and isinstance(v2, np.ndarray):
            both_nan = np.isnan(v1) & np.isnan(v2)
            diff = ~both_nan & (v1 != v2)
            return "SAME" if np.count_nonzero(diff) == 0 else "DIFFERENT"

        if isinstance(v1, float) and isinstance(v2, float):
            if math.isnan(v1) and math.isnan(v2):
                return "SAME"

        return "SAME" if v1 == v2 else "DIFFERENT"

    # Read files
    d1 = _read_bin_to_dict(file1)
    d2 = _read_bin_to_dict(file2)

    # Order fields
    keys = sorted(set(d1) | set(d2))
    rows = []

    # Loop over fields [rows]
    for k in keys:
        r1 = d1.get(k)
        r2 = d2.get(k)

        if r1 is None or r2 is None:
            rows.append((k, "MISSING", _summarize(r1), _summarize(r2)))
        else:
            status = _compare(r1, r2)
            rows.append((k, status, _summarize(r1), _summarize(r2)))

    # Construct output
    headers = ("Field", "Status", f"File 1 ({mode})", f"File 2 ({mode})")
    colw = [
        max(len(h), *(len(str(r[i])) for r in rows))
        for i, h in enumerate(headers)
    ]

    lines = []
    lines.append(f"Comparing ({mode}):\n  File 1: {file1}\n  File 2: {file2}\n")
    lines.append(
        f"{headers[0]:<{colw[0]}}  {headers[1]:<{colw[1]}}  "
        f"{headers[2]:<{colw[2]}}  {headers[3]:<{colw[3]}}"
    )
    lines.append("-" * (sum(colw) + 6))

    for r in rows:
        lines.append(
            f"{r[0]:<{colw[0]}}  {r[1]:<{colw[1]}}  "
            f"{r[2]:<{colw[2]}}  {r[3]:<{colw[3]}}"
        )

    # Print output (to file or terminal)
    if out_file:
        with open(out_file, "w") as f:
            f.write("\n".join(lines))
        print(f'Output written to file: {out_file}.')
    else:
        print("\n".join(lines))