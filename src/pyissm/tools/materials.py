"""
Material property functions for ISSM

This module contains functions to compute ice material properties.
"""

import numpy as np
import warnings


def paterson(temperature):

    """
    Compute ice rigidity using the Paterson parameterization.

    Parameters
    ----------
    temperature : array_like
        Temperature value(s) in Kelvin. Must be non-negative. Scalar or array-like
        inputs are accepted and will be converted to a NumPy array.

    Returns
    -------
    rigidity : ndarray
        Array of computed rigidity values with the same shape as ``temperature``.
        The return dtype is float. Values that are computed as negative are
        replaced by a floor value of 1e6.

    Raises
    ------
    RuntimeError
        If any element of ``temperature`` is negative (temperatures must be given
        in Kelvin).

    Notes
    -----
    The implementation converts input temperatures from Kelvin to Celsius and
    then evaluates a piecewise cubic polynomial fit over several temperature
    intervals. A temperature-dependent shift is applied within each interval
    before polynomial evaluation. A global scale factor (1e8) is multiplied by
    the polynomial result. To ensure physical plausibility, any negative
    computed rigidity values are replaced with 1e6.

    The function preserves the shape of the input (broadcasting rules of NumPy
    apply) and returns an array of the same shape.

    Warnings
    --------
    DeprecationWarning
        This function is outdated. Use :func:`cuffey` instead.
    """

    # Deprecation warning
    warnings.warn('pyissm.tools.materials.paterson: paterson() is outdated. Use pyissm.tools.materials.cuffey instead.', DeprecationWarning)
    
    # Ensure input is numpy array
    temperature = np.asarray(temperature)

    # Error checks
    if np.any(temperature < 0.0):
        raise RuntimeError("input temperature should be in Kelvin (positive)")
    
    # Convert to Celsius
    T = temperature - 273.15

    # Initialize rigidity array
    rigidity = np.zeros_like(T, dtype=float)

    pieces = [
        (-np.inf, -45,  50, [-0.000292866376675, 0.011672640664130, -0.325004442485481, 6.524779401948101]),
        (-45,     -40,  45, [-0.000292866376675, 0.007279645014004, -0.230243014094813, 5.154964909039554]),
        (-40,     -35,  40, [0.000072737147457, 0.002886649363879, -0.179411542205399, 4.149132666831214]),
        (-35,     -30,  35, [-0.000086144770023, 0.003977706575736, -0.145089762507325, 3.333333333333331]),
        (-30,     -25,  30, [-0.000043984685769, 0.002685535025386, -0.111773554501713, 2.696559088937191]),
        (-25,     -20,  25, [-0.000029799523463, 0.002025764738854, -0.088217055680511, 2.199331606342181]),
        (-20,     -15,  20, [0.000136920904777, 0.001578771886910, -0.070194372551690, 1.805165505978111]),
        (-15,     -10,  15, [-0.000899763781026, 0.003632585458564, -0.044137585824322, 1.510778053489523]),
        (-10,      -5,  10, [0.001676964325070, -0.009863871256831, -0.075294014815659, 1.268434288203714]),
        ( -5,      -2,   5, [-0.003748937622487, 0.015290593619213, -0.048160403003748, 0.854987973338348]),
        ( -2,   np.inf,  2, [-0.003748937622488, -0.018449844983174, -0.057638157095631, 0.746900791092860]),
    ]

    # Evaluate each piece
    scale = 1e8
    for Tmin, Tmax, shift, (c3, c2, c1, c0) in pieces:
        mask = (T >= Tmin) & (T < Tmax)
        if np.any(mask):
            x = T[mask] + shift
            rigidity[mask] = scale * (c3*x**3 + c2*x**2 + c1*x + c0)

    # Enforce positivity
    rigidity = np.where(rigidity < 0, 1e6, rigidity)

    return rigidity


def cuffey(temperature):
    """
    Compute ice rigidity using the Cuffey parameterization.

    Parameters
    ----------
    temperature : array_like
        Temperature value(s) in Kelvin. Must be non-negative. Scalar or array-like
        inputs are accepted and will be converted to a NumPy array.

    Returns
    -------
    rigidity : ndarray
        Array of computed rigidity values with the same shape as ``temperature``.
        The return dtype is float. Values that are computed as negative are
        replaced by a floor value of 1e6.

    Raises
    ------
    RuntimeError
        If any element of ``temperature`` is negative (temperatures must be given
        in Kelvin).

    Notes
    -----
    The implementation converts input temperatures from Kelvin to Celsius and
    then evaluates a piecewise cubic polynomial fit over several temperature
    intervals. A temperature-dependent shift is applied within each interval
    before polynomial evaluation. A global scale factor (1e8) is multiplied by
    the polynomial result. To ensure physical plausibility, any negative
    computed rigidity values are replaced with 1e6.

    The function preserves the shape of the input (broadcasting rules of NumPy
    apply) and returns an array of the same shape.
    """

    # Ensure input is numpy array
    temperature = np.asarray(temperature)

    # Error checks
    if np.any(temperature < 0.0):
        raise RuntimeError("input temperature should be in Kelvin (positive)")

    # Convert to Celsius
    T = temperature - 273.15

    # Initialize rigidity array
    rigidity = np.zeros_like(T, dtype=float)

    # Define piecewise ranges and polynomial coefficients
    # Each entry: (T_min, T_max, shift, [c3, c2, c1, c0])
    pieces = [
        (-np.inf, -45, 50, [-0.000396645116301,  0.013345579471334, -0.356868703259105, 7.272363035371383]),
        (-45, -40, 45,     [-0.000396645116301,  0.007395902726819, -0.253161292268336, 5.772078366321591]),
        (-40, -35, 40,     [ 0.000408322072669,  0.001446225982305, -0.208950648722716, 4.641588833612773]),
        (-35, -30, 35,     [-0.000423888728124,  0.007571057072334, -0.163864233449525, 3.684031498640382]),
        (-30, -25, 30,     [ 0.000147154327025,  0.001212726150476, -0.119945317335478, 3.001000667185614]),
        (-25, -20, 25,     [-0.000193435838672,  0.003420041055847, -0.096781481303861, 2.449986525148220]),
        (-20, -15, 20,     [ 0.000219771255067,  0.000518503475772, -0.077088758645767, 2.027400665191131]),
        (-15, -10, 15,     [-0.000653438900191,  0.003815072301777, -0.055420879758021, 1.682390865739973]),
        (-10, -5, 10,      [ 0.000692439419762, -0.005986511201093, -0.066278074254598, 1.418983411970382]),
        (-5, -2,  5,       [-0.000132282004110,  0.004400080095332, -0.074210229783403, 1.024485188140279]),
        (-2, np.inf, 2,    [-0.000132282004110,  0.003209542058346, -0.051381363322371, 0.837883605537096])
    ]

    # Evaluate each piece
    scale = 1e8
    for Tmin, Tmax, shift, (c3, c2, c1, c0) in pieces:
        mask = (T >= Tmin) & (T < Tmax)
        if np.any(mask):
            x = T[mask] + shift
            rigidity[mask] = scale * (c3*x**3 + c2*x**2 + c1*x + c0)

    # Enforce positivity
    rigidity = np.where(rigidity < 0, 1e6, rigidity)

    return rigidity
