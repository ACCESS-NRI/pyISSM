"""
Utility functions for ISSM

This module contains various utility functions that are used throughout the ISSM codebase.
"""

import numpy as np

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
        raise ValueError(f"Invalid input_units '{input_units}'. Must be one of {accepted_units}")

    if output_units not in accepted_units:
        raise ValueError(f"Invalid output_units '{output_units}'. Must be one of {accepted_units}")

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
        raise ValueError("""Requested unit conversion currently not supported. Available conversions include:
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
