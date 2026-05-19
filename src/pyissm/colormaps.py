"""
Custom colormaps for ISSM visualizations.

Ported from ISSM src/m/plot/colormaps/
"""

import numpy as np
import matplotlib as mpl
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, hsv_to_rgb


def seacolor(n=256):
    """Sea colormap (deep navy to pale blue).

    Parameters
    ----------
    n : int
        Number of color levels. Default 256.

    Returns
    -------
    np.ndarray, shape (3, n)
        RGB color array.
    """
    J = [
        0.0392, 0, 0.4745,
        0.1020, 0, 0.5373,
        0.1020, 0, 0.5373,
        0.1490, 0, 0.5961,
        0.1490, 0, 0.5961,
        0.1059, 0.0118, 0.6510,
        0.1059, 0.0118, 0.6510,
        0.0627, 0.0235, 0.7059,
        0.0627, 0.0235, 0.7059,
        0.0196, 0.0353, 0.7569,
        0.0196, 0.0353, 0.7569,
        0, 0.0549, 0.7961,
        0, 0.0549, 0.7961,
        0, 0.0863, 0.8235,
        0, 0.0863, 0.8235,
        0, 0.1176, 0.8471,
        0, 0.1176, 0.8471,
        0, 0.1529, 0.8745,
        0, 0.1529, 0.8745,
        0.0471, 0.2667, 0.9059,
        0.0471, 0.2667, 0.9059,
        0.1020, 0.4000, 0.9412,
        0.1020, 0.4000, 0.9412,
        0.0745, 0.4588, 0.9569,
        0.0745, 0.4588, 0.9569,
        0.0549, 0.5216, 0.9765,
        0.0549, 0.5216, 0.9765,
        0.0824, 0.6196, 0.9882,
        0.0824, 0.6196, 0.9882,
        0.1176, 0.6980, 1.0000,
        0.1176, 0.6980, 1.0000,
        0.1686, 0.7294, 1.0000,
        0.1686, 0.7294, 1.0000,
        0.2157, 0.7569, 1.0000,
        0.2157, 0.7569, 1.0000,
        0.2549, 0.7843, 1.0000,
        0.2549, 0.7843, 1.0000,
        0.3098, 0.8235, 1.0000,
        0.3098, 0.8235, 1.0000,
        0.3686, 0.8745, 1.0000,
        0.3686, 0.8745, 1.0000,
        0.5412, 0.8902, 1.0000,
        0.5412, 0.8902, 1.0000,
        0.7373, 0.9020, 1.0000,
    ]
    length = int(len(J) / 3)
    J = np.array(J).reshape(length, 3)
    a = np.linspace(1, length, n)
    b = np.arange(1, length + 1)
    return np.array([np.interp(a, b, J[:, i]) for i in range(3)])


def landcolor(n=256):
    """Land colormap (green to brown).

    Parameters
    ----------
    n : int
        Number of color levels. Default 256.

    Returns
    -------
    np.ndarray, shape (3, n)
        RGB color array.
    """
    J = [
        0.095678, 0.53427, 0.21682,
        0.15785, 0.5979, 0.23274,
        0.21286, 0.64673, 0.2514,
        0.26411, 0.68789, 0.27268,
        0.32959, 0.72416, 0.31308,
        0.39794, 0.75695, 0.36038,
        0.46153, 0.7871, 0.40624,
        0.52108, 0.81516, 0.45135,
        0.57702, 0.84152, 0.49547,
        0.62973, 0.86645, 0.53891,
        0.67946, 0.89016, 0.58187,
        0.72647, 0.91282, 0.62427,
        0.77095, 0.93455, 0.66619,
        0.81306, 0.95546, 0.70772,
        0.85292, 0.97563, 0.7489,
        0.89066, 0.99514, 0.78976,
        0.88379, 0.98595, 0.77038,
        0.86389, 0.96758, 0.73236,
        0.84615, 0.94972, 0.69623,
        0.8303, 0.93233, 0.66186,
        0.81612, 0.91536, 0.6291,
        0.80341, 0.8988, 0.59784,
        0.79201, 0.8826, 0.56795,
        0.78191, 0.86676, 0.53946,
        0.7729, 0.85123, 0.51224,
        0.76479, 0.83602, 0.48615,
        0.75747, 0.8211, 0.46111,
        0.75084, 0.80645, 0.43704,
        0.74506, 0.79206, 0.41414,
        0.73981, 0.77792, 0.39211,
        0.73501, 0.76401, 0.37089,
        0.73068, 0.75033, 0.35052,
        0.72683, 0.73685, 0.33106,
        0.72042, 0.72074, 0.31228,
        0.71032, 0.70085, 0.29417,
        0.69761, 0.67821, 0.27694,
        0.68489, 0.65558, 0.26026,
        0.67235, 0.63313, 0.24418,
        0.65997, 0.61082, 0.22889,
        0.64775, 0.58874, 0.21406,
        0.63568, 0.56689, 0.19983,
        0.62376, 0.54527, 0.18622,
        0.61197, 0.52391, 0.17299,
        0.60033, 0.50283, 0.16046,
        0.58881, 0.48203, 0.14832,
        0.57742, 0.46151, 0.13667,
        0.56616, 0.44133, 0.12555,
        0.55502, 0.4214, 0.11472,
        0.54398, 0.4019, 0.10456,
        0.53306, 0.38266, 0.094633,
        0.52226, 0.36382, 0.085242,
        0.51155, 0.3453, 0.076179,
        0.50095, 0.32714, 0.067515,
        0.49045, 0.30938, 0.059259,
        0.48005, 0.29193, 0.051294,
        0.46973, 0.27495, 0.043796,
        0.45951, 0.25823, 0.0365,
        0.44938, 0.24206, 0.029715,
        0.43934, 0.22609, 0.023063,
        0.42938, 0.21074, 0.016949,
        0.41951, 0.19556, 0.010917,
        0.40971, 0.18105, 0.0054326,
        0.4, 0.16667, 0,
    ]
    length = int(len(J) / 3)
    J = np.array(J).reshape(length, 3)
    a = np.linspace(1, length, n)
    b = np.arange(1, length + 1)
    return np.array([np.interp(a, b, J[:, i]) for i in range(3)])


def ibcao(nsea, nland):
    """IBCAO colormap for polar regions.

    Parameters
    ----------
    nsea : int
        Number of sea colors.
    nland : int
        Number of land colors.

    Returns
    -------
    np.ndarray, shape (nsea + nland, 3)
        RGB color array.
    """
    Jsea = [
        0.18039, 0.29020, 0.57255,
        0.18039, 0.29020, 0.57255,
        0.05882, 0.44314, 0.65490,
        0.05882, 0.44314, 0.65490,
        0.02745, 0.49804, 0.73725,
        0.02745, 0.49804, 0.73725,
        0.01176, 0.54510, 0.78824,
        0.01176, 0.54510, 0.78824,
        0.00784, 0.63529, 0.83922,
        0.00784, 0.63529, 0.83922,
        0.06667, 0.71765, 0.86667,
        0.06667, 0.71765, 0.86667,
        0.17647, 0.75294, 1.00000,
        0.17647, 0.75294, 1.00000,
        0.23529, 0.76471, 0.85882,
        0.23529, 0.76471, 0.85882,
        0.24314, 0.76471, 0.83922,
        0.24314, 0.76471, 0.83922,
        0.25882, 0.76078, 0.81176,
        0.25882, 0.76078, 0.81176,
        0.27451, 0.76078, 0.76078,
        0.27451, 0.76078, 0.76078,
        0.41961, 0.78431, 0.74902,
        0.41961, 0.78431, 0.74902,
        0.60000, 0.83137, 0.74902,
        0.60000, 0.83137, 0.74902,
    ]
    Jland = [
        0.85098, 0.84314, 0.30588,
        0.85098, 0.84314, 0.30588,
        0.93333, 0.89020, 0.41961,
        0.93333, 0.89020, 0.41961,
        0.93725, 0.80784, 0.35686,
        0.93725, 0.80784, 0.35686,
        0.89804, 0.74510, 0.31765,
        0.89804, 0.74510, 0.31765,
        0.85098, 0.63922, 0.21961,
        0.85098, 0.63922, 0.21961,
        0.75686, 0.55294, 0.22353,
        0.75686, 0.55294, 0.22353,
        0.71765, 0.50980, 0.22353,
        0.71765, 0.50980, 0.22353,
        0.68627, 0.48235, 0.21961,
        0.68627, 0.48235, 0.21961,
        0.65490, 0.45882, 0.21569,
        0.65490, 0.45882, 0.21569,
        0.58824, 0.39608, 0.20392,
        0.58824, 0.39608, 0.20392,
        1.00000, 1.00000, 1.00000,
    ]
    lsea = int(len(Jsea) / 3)
    Jsea = np.array(Jsea).reshape(lsea, 3)
    a = np.linspace(1, lsea, nsea)
    b = np.arange(1, lsea + 1)
    ysea = np.array([np.interp(a, b, Jsea[:, i]) for i in range(3)])

    lland = int(len(Jland) / 3)
    Jland = np.array(Jland).reshape(lland, 3)
    a = np.linspace(1, lland, nland)
    b = np.arange(1, lland + 1)
    yland = np.array([np.interp(a, b, Jland[:, i]) for i in range(3)])

    return np.concatenate((ysea, yland), axis=1).T


def demmap(ncolors, minZ, maxZ, colorscheme='dem'):
    """DEM colormap concatenating sea and land colors.

    Parameters
    ----------
    ncolors : int
        Total number of colors.
    minZ : float
        Minimum elevation (negative for seafloor).
    maxZ : float
        Maximum elevation.
    colorscheme : str
        'dem' (default) or 'ibcao'.

    Returns
    -------
    matplotlib.colors.ListedColormap

    Examples
    --------
    cmap = demmap(50, -300, 1200)
    cmap = demmap(50, -300, 1200, 'ibcao')
    """
    if minZ == maxZ:
        maxZ = minZ + 1

    if minZ >= 0:
        nsea, nland = 0, ncolors
    elif maxZ <= 0:
        nsea, nland = ncolors, 0
    else:
        maxminratio = maxZ / abs(minZ)
        n1 = int(np.floor(ncolors / 2))
        n2 = int(np.ceil(ncolors / 2))
        if maxminratio > 1:
            sea_arr = np.arange(1, n1 + 1)
            land_arr = np.arange(ncolors - 1, n2 - 1, -1)
        else:
            land_arr = np.arange(1, n1 + 1)
            sea_arr = np.arange(ncolors - 1, n2 - 1, -1)
        ratio = land_arr / sea_arr
        idx = int(np.argmin(np.abs(ratio - maxminratio) / maxminratio))
        nsea = int(sea_arr[idx])
        nland = int(land_arr[idx])

    if colorscheme.lower() == 'dem':
        colors = np.concatenate((seacolor(nsea), landcolor(nland) ** 1.3), axis=1).T
    elif colorscheme.lower() == 'ibcao':
        colors = ibcao(nsea, nland)
    else:
        raise ValueError(f"demmap: unknown colorscheme '{colorscheme}', use 'dem' or 'ibcao'")

    return ListedColormap(colors)


def get_colormap(name, alpha=1):
    """Get a matplotlib colormap by name, including ISSM custom maps.

    Supported custom names:
        'Rignot'         — modified HSV rainbow scaled by velocity
        'demmap(n,z0,z1)'  — calls demmap(n, z0, z1)
        Any standard matplotlib colormap name

    Parameters
    ----------
    name : str or matplotlib.colors.Colormap
        Colormap name or existing colormap object.
    alpha : float
        Exponent used for 'Rignot' colormap. Default 1.

    Returns
    -------
    matplotlib.colors.Colormap
    """
    if isinstance(name, mpl.colors.Colormap):
        return name
    if name in mpl.colormaps:
        return mpl.colormaps[name]
    if not name:
        return mpl.colormaps['viridis']
    if name == 'Rignot':
        c = np.array((np.linspace(0, 1, 128, False), np.ones(128), np.ones(128))).T
        c[:, 1] = np.clip(0.1 + c[:, 0] ** (1 / alpha), 0, 1)
        return ListedColormap(hsv_to_rgb(c))
    try:
        return ListedColormap(eval(name))
    except Exception:
        raise ValueError(
            f"get_colormap: '{name}' is not a recognized colormap. "
            "Use a matplotlib name, 'Rignot', or a demmap/ibcao expression."
        )


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    """Truncate a colormap within normalized limits [minval, maxval].

    Parameters
    ----------
    cmap : str or matplotlib.colors.Colormap
    minval : float
        Lower bound, normalized [0, 1].
    maxval : float
        Upper bound, normalized [0, 1].
    n : int
        Number of levels in the new colormap.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap

    Example
    -------
    cmap = truncate_colormap('viridis', 0.2, 0.8)
    """
    if isinstance(cmap, str):
        cmap = mpl.colormaps[cmap]
    return LinearSegmentedColormap.from_list(
        f'trunc({cmap.name},{minval:.2f},{maxval:.2f})',
        cmap(np.linspace(minval, maxval, n))
    )
