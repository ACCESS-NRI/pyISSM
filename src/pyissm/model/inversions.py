"""
Inversion utilities for ISSM models.

Provides the canonical lists of supported control parameters and cost
functions, plus a helper for converting numeric cost-function IDs to
their string names (as expected by ISSM's marshalling layer).
"""

# Maps numeric cost-function IDs to the string names used internally by ISSM.
# Velocity misfits (101–105), thickness misfit (201),
# gradient regularisation terms (501–505).
COST_FUNCTION_NAMES = {
    101: 'SurfaceAbsVelMisfit',
    102: 'SurfaceRelVelMisfit',
    103: 'SurfaceLogVelMisfit',
    104: 'SurfaceLogVxVyMisfit',
    105: 'SurfaceAverageVelMisfit',
    201: 'ThicknessAbsMisfit',
    501: 'DragCoefficientAbsGradient',
    502: 'RheologyBbarAbsGradient',
    503: 'ThicknessAbsGradient',
    504: 'ThicknessAlongGradient',
    505: 'ThicknessAcrossGradient',
}


def supported_controls():
    """Return the list of field names that can be used as inversion controls.

    Returns
    -------
    list of str
        Control parameter names accepted by md.inversion.control_parameters.

    Examples
    --------
    >>> from pyissm.model import inversions
    >>> inversions.supported_controls()
    ['BalancethicknessThickeningRate', 'FrictionCoefficient', ...]
    """
    return [
        'BalancethicknessThickeningRate',
        'FrictionCoefficient',
        'FrictionC',
        'FrictionAs',
        'MaterialsRheologyBbar',
        'DamageDbar',
        'Vx',
        'Vy',
        'Thickness',
        'BalancethicknessOmega',
        'BalancethicknessApparentMassbalance',
        'MaterialsRheologyB',
    ]


def supported_cost_functions():
    """Return the list of numeric IDs for supported inversion cost functions.

    Returns
    -------
    list of int
        Cost-function IDs accepted by md.inversion.cost_functions.
        Use :func:`cost_function_name` to convert an ID to its string name.

    Examples
    --------
    >>> from pyissm.model import inversions
    >>> inversions.supported_cost_functions()
    [101, 102, 103, 104, 105, 201, 501, 502, 503, 504, 505]
    """
    return list(COST_FUNCTION_NAMES.keys())


def cost_function_name(cost_function_id):
    """Convert a numeric cost-function ID to its ISSM string name.

    Parameters
    ----------
    cost_function_id : int or list of int
        One or more cost-function IDs (e.g. 101, [101, 503]).

    Returns
    -------
    str or list of str
        Corresponding ISSM string name(s).

    Raises
    ------
    KeyError
        If an ID is not in :data:`COST_FUNCTION_NAMES`.

    Examples
    --------
    >>> inversions.cost_function_name(101)
    'SurfaceAbsVelMisfit'
    >>> inversions.cost_function_name([101, 503])
    ['SurfaceAbsVelMisfit', 'ThicknessAbsGradient']
    """
    if isinstance(cost_function_id, int):
        return COST_FUNCTION_NAMES[cost_function_id]
    return [COST_FUNCTION_NAMES[cf] for cf in cost_function_id]
