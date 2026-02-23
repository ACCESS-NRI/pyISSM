"""
Various utility functions for ISSM model classes.
"""

import numpy as np
import os
import re
from operator import attrgetter
from pyissm.model.classes import basalforcings, calving, friction, smb, frontalforcings

"""
Functions for formatting and displaying model fields in ISSM classes. Taken from $ISSM_DIR/src/m/miscellaneous/fielddisplay.py.
"""

def fielddisplay(md, name, comment):
    """
    Display a model field with formatted output for ISSM class representations.

    This function extracts a field from a model object and formats it for display
    in the __repr__ method of ISSM classes. It handles various data types and
    provides consistent formatting across all model components.

    Parameters
    ----------
    md : :class:`object`
        The model object containing the field to display.
    name : :class:`str`
        The name of the field attribute to extract and display.
    comment : :class:`str`
        Descriptive comment explaining the field's purpose and units.

    Returns
    -------
    :class:`str`
        Formatted string representation of the field suitable for display.

    Examples
    --------
    .. code-block:: python

        >>> fielddisplay(md, 'thickness', 'ice thickness [m]')
        '         thickness        : (100, 200)    -- ice thickness [m]'
        
        >>> fielddisplay(md, 'verbose', 'verbose output flag')
        '         verbose          : True          -- verbose output flag'
    """

    #get field
    field = getattr(md, name)

    #disp corresponding line as a function of field type (offset set as 9 spaces)
    return parsedisplay("         ", name, field, comment)


def parsedisplay(offset, name, field, comment):
    """
    Parse and format a field for display based on its data type.

    This function determines the appropriate formatting for a field based on its
    data type (string, numeric, array, boolean, dictionary, list, etc.) and
    returns a consistently formatted string representation.

    Parameters
    ----------
    offset : :class:`str`
        String of spaces used for indentation in the output.
    name : :class:`str`
        The name of the field being displayed.
    field : any
        The field value to be formatted for display.
    comment : :class:`str`
        Descriptive comment explaining the field's purpose and units.

    Returns
    -------
    :class:`str`
        Formatted string representation of the field.

    Notes
    -----
    Handles the following data types:

        - :class:`str`: Shows string value (truncated if > 30 characters)
        - :class:`int`, :class:`float`: Shows numeric value
        - :class:`numpy.ndarray`: Shows array shape
        - :class:`bool`: Shows True/False
        - :class:`dict`: Shows nested dictionary structure
        - :class:`list`, :class:`tuple`: Shows container contents (up to 5 elements)
        - ``None`: Shows "None"
        - Other types: Shows "not displayed"

    Examples
    --------
    .. code-block:: python

            >>> parsedisplay('   ', 'temperature', 273.15, 'temperature [K]')
        '   temperature        : 273.15        -- temperature [K]'
    """
    #string
    if isinstance(field, str):
        if len(field) > 30:
            string = displayunit(offset, name, "not displayed", comment)
        else:
            string = displayunit(offset, name, "'%s'" % field, comment)

    #numeric
    elif isinstance(field, (int, float)):
        string = displayunit(offset, name, str(field), comment)

    #matrix
    elif isinstance(field, np.ndarray):
        string = displayunit(offset, name, str(field.shape), comment)

    #logical
    elif isinstance(field, bool):
        if field:
            string = displayunit(offset, name, "True", comment)
        else:
            string = displayunit(offset, name, "False", comment)

    #dictionary
    elif isinstance(field, dict):
        string = dict_display(offset, name, field, comment)

    #list or tuple
    elif isinstance(field, (list, tuple)):
        string = list_display(offset, name, field, comment)

    #None
    elif field is None:
        string = displayunit(offset, name, "None", comment)

    else:
        string = displayunit(offset, name, "not displayed", comment)

    return string


def dict_display(offset, name, field, comment):
    """
    Format a dictionary field for hierarchical display.

    This function formats dictionary fields by displaying the main dictionary
    indicator and then recursively formatting each key-value pair with
    increased indentation to show the hierarchical structure.

    Parameters
    ----------
    offset : :class:`str`
        String of spaces used for base indentation in the output.
    name : :class:`str`
        The name of the dictionary field being displayed.
    field : dict
        The dictionary to be formatted for display.
    comment : :class:`str`
        Descriptive comment explaining the dictionary's purpose.

    Returns
    -------
    :class:`str`
        Formatted multi-line string representation of the dictionary structure.

    Notes
    -----
    Empty dictionaries are displayed as 'N/A'. Non-empty dictionaries show
    '{dictionary}' on the first line followed by indented key-value pairs.

    Examples
    --------
    .. code-block:: python

        >>> dict_display('   ', 'params', {'a': 1, 'b': 2}, 'parameters')
        '   params             : {dictionary}   -- parameters
            a               : 1
            b               : 2'
    """
    if field:
        string = displayunit(offset, name, '{dictionary}', comment) + '\n'
        offset += '   '

        for structure_field, sfield in list(field.items()):
            string += parsedisplay(offset, str(structure_field), sfield, '') + '\n'

        if string and string[-1] == '\n':
            string = string[:-1]

    else:
        string = displayunit(offset, name, 'N/A', comment)

    return string


def list_display(offset, name, field, comment):
    """
    Format a list or tuple field for compact display.

    This function formats list and tuple fields by showing their contents
    in a compact format. For small containers (< 5 elements), it shows all
    elements. For larger containers, it shows the size.

    Parameters
    ----------
    offset : :class:`str`
        String of spaces used for indentation in the output.
    name : :class:`str`
        The name of the list/tuple field being displayed.
    field : list or tuple
        The list or tuple to be formatted for display.
    comment : :class:`str`
        Descriptive comment explaining the container's purpose.

    Returns
    -------
    :class:`str`
        Formatted string representation of the list or tuple.

    Notes
    -----
        
    - Small containers show individual elements: [1, 2, 3] or ('a', 'b', 'c').
    - Large containers show size: [100x1] or (50x1).
    - Only simple data types (:class:`str`, :class:`bool`, :class:`int`, :class:`float`) are shown individually.

    Examples
    --------
    .. code-block:: python

        >>> list_display('   ', 'coords', [1.0, 2.0, 3.0], 'coordinates')
        '   coords             : [1.0, 2.0, 3.0] -- coordinates'
        
        >>> list_display('   ', 'data', list(range(100)), 'large dataset')
        '   data               : [100x1]         -- large dataset'
    """

    #initialization
    if isinstance(field, list):
        sbeg = '['
        send = ']'
    elif isinstance(field, tuple):
        sbeg = '('
        send = ')'
    string = sbeg

    #go through the cell and fill string
    if len(field) < 5:
        for fieldi in field:
            if isinstance(fieldi, str):
                string += "'%s', " % fieldi
            elif isinstance(fieldi, (bool, int, float)):
                string += "%s, " % str(fieldi)
            else:
                string = sbeg
                break

    if string == sbeg:
        string = "%s%dx1%s" % (sbeg, len(field), send)
    else:
        string = string[:-1] + send

    #call displayunit
    return displayunit(offset, name, string, comment)


def displayunit(offset, name, characterization, comment):
    """
    Format a single field line with consistent spacing and alignment.

    This is the core formatting function that creates the final display line
    with proper spacing, alignment, and comment formatting. It handles name
    truncation, characterization formatting, and multi-line comments.

    Parameters
    ----------
    offset : :class:`str`
        String of spaces used for indentation.
    name : :class:`str`
        The field name (truncated if > 23 characters).
    characterization : :class:`str`
        The string representation of the field value.
    comment : :class:`str` or list of str
        Descriptive comment(s) explaining the field. Can be a single string
        or a list of strings for multi-line comments.

    Returns
    -------
    :class:`str`
        Formatted display line with proper spacing and alignment.

    Notes
    -----
    Format: "{offset}{name:23}: {characterization:15} -- {comment}"
    
    Special handling:

        - Names > 23 chars are truncated with "..."
        - Characterizations > 15 chars are truncated with "..."
        - Empty/NaN values are shown as "N/A"
        - List comments create multi-line output

    Examples
    --------
    .. code-block:: python

        >>> displayunit('   ', 'temperature', '273.15', 'temperature [K]')
        '   temperature        : 273.15         -- temperature [K]'
        
        >>> displayunit('   ', 'very_long_field_name', '42', 'description')
        '   very_long_field_n...: 42             -- description'
    """

    #take care of name
    if len(name) > 23:
        name = "%s..." % name[:20]

    #take care of characterization
    if characterization in ["''", '""', 'nan', np.nan, 'NaN', "[0x1]"]:
        characterization = "N/A"

    if len(characterization) > 15:
        characterization = "%s..." % characterization[:12]

    #print
    if not comment:
        string = "%s% - 23s: % - 15s" % (offset, name, characterization)
    else:
        if isinstance(comment, str):
            string = "%s% - 23s: % - 15s -- %s" % (offset, name, characterization, comment)
        elif isinstance(comment, list):
            string = "%s% - 23s: % - 15s -- %s" % (offset, name, characterization, comment[0])
            for commenti in comment:
                string += "\n%s% - 23s  % - 15s    %s" % (offset, '', '', commenti)
        else:
            raise RuntimeError("fielddisplay error message: format for comment not supported yet")

    return string

def getlongestfieldname(self):
    """
    Find the longest field name in an object's attributes.

    This utility function iterates through all attributes of an object
    and returns the length of the longest field name. This can be used
    for dynamic formatting and alignment purposes.

    Parameters
    ----------
    self : object
        The object whose attribute names will be examined.

    Returns
    -------
    :class:`int`
        The length of the longest field name in the object's __dict__.

    Notes
    -----
    Only considers attributes in the object's __dict__, not inherited
    attributes or properties.

    Examples
    --------
    ..code-block:: python

        >>> class Example:
        ...     def __init__(self):
        ...         self.a = 1
        ...         self.very_long_field_name = 2
        ...         self.b = 3
        >>> obj = Example()
        >>> getlongestfieldname(obj)
        20
    
    """

    maxlength = 0
    for key in self.__dict__.keys():
        length = len(key)
        if length > maxlength:
            maxlength = length

    return maxlength

def marshall_inversion_cost_functions(cost_functions):
    """
    Map inversion cost function codes to their corresponding string names.
    
    This function converts integer cost function codes used in ISSM inversion
    routines to their human-readable string representations. It supports both
    single cost function codes and lists of multiple codes.
    
    Parameters
    ----------
    cost_functions : :class:`int` or list of :class:`int`
        Cost function code(s) to be mapped. Supported codes include:

        - 101: SurfaceAbsVelMisfit
        - 102: SurfaceRelVelMisfit  
        - 103: SurfaceLogVelMisfit
        - 104: SurfaceLogVxVyMisfit
        - 105: SurfaceAverageVelMisfit
        - 201: ThicknessAbsMisfit
        - 501: DragCoefficientAbsGradient
        - 502: RheologyBbarAbsGradient
        - 503: ThicknessAbsGradient
        - 504: ThicknessAlongGradient
        - 505: ThicknessAcrossGradient
    
    Returns
    -------
    list of :class:`str`
        List containing the string name(s) corresponding to the input cost
        function code(s).
    
    Examples
    --------
    .. code-block:: python

        >>> map_inversion_cost_functions(101)
        ['SurfaceAbsVelMisfit']
        
        >>> map_inversion_cost_functions([101, 201, 501])
        ['SurfaceAbsVelMisfit', 'ThicknessAbsMisfit', 'DragCoefficientAbsGradient']
    """

    ## Define dictionary of cost functions
    cfDict = {101: 'SurfaceAbsVelMisfit',
            102: 'SurfaceRelVelMisfit',
            103: 'SurfaceLogVelMisfit',
            104: 'SurfaceLogVxVyMisfit',
            105: 'SurfaceAverageVelMisfit',
            201: 'ThicknessAbsMisfit',
            501: 'DragCoefficientAbsGradient',
            502: 'RheologyBbarAbsGradient',
            503: 'ThicknessAbsGradient',
            504: 'ThicknessAlongGradient',
            505: 'ThicknessAcrossGradient'}

    ## Marshall cost functions
    if isinstance(cost_functions, int):
        return [cfDict[cost_functions]]
    else:
        return [cfDict[cf] for cf in cost_functions]

def supported_inversion_control_parameters():
    """
    Return a list of supported inversion control parameters.

    Returns
    -------
    list of :class:`str`
        List of supported inversion control parameter names.
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
        'MaterialsRheologyB'
    ]

def supported_inversion_cost_functions():
    """
    Return a list of supported inversion cost function codes.

    Returns
    -------
    list of :class:`int`
        List of supported inversion cost function codes.
    """
    return list(range(101, 105 + 1)) + [201] + list(range(501, 508 + 1)) + [510] + list(range(601, 604 + 1))

def supported_analyses():
    """
    Return a list of supported analyses.

    Returns
    -------
    list of :class:`str`
        List of supported analyses.
    """
    return [
        'DefaultAnalysis',
        'RecoveryAnalysis',
        'StressbalanceAnalysis',
        'StressbalanceVerticalAnalysis',
        'GLheightadvectionAnalysis',
        'MasstransportAnalysis',
        'ThermalAnalysis',
        'EnthalpyAnalysis',
        'AdjointBalancethicknessAnalysis',
        'BalancethicknessAnalysis',
        'Balancethickness2Analysis',
        'BalancethicknessSoftAnalysis',
        'BalancevelocityAnalysis',
        'DamageEvolutionAnalysis',
        'LoveAnalysis',
        'EsaAnalysis',
        'SealevelchangeAnalysis',
        'FreeSurfaceBaseAnalysis',
        'FreeSurfaceTopAnalysis',
        'LevelsetAnalysis',
        'DebrisAnalysis',
        'L2ProjectionBaseAnalysis',
        'ExtrudeFromBaseAnalysis',
        'ExtrudeFromTopAnalysis'
    ]

def supported_stochastic_forcings(return_dict = False):
    """
    Return supported stochastic forcings.
    
    Parameters
    ----------
    return_dict : :class:`bool`, default=False
        If True, return the full mapping dict; 
        otherwise, return just the list of keys.

    Returns
    -------
    :class:`dict` or list of :class:`str`
        Mapping of supported stochastic forcing fields to their
        corresponding parameter modules, or list of field names.

    Notes
    -----
    Maintain legacy naming for compatibility with MATLAB version
    """

    structure = {
        'BasalforcingsDeepwaterMeltingRatearma': basalforcings.lineararma,
        'BasalforcingsSpatialDeepwaterMeltingRate': basalforcings.spatiallinear,
        'DefaultCalving': calving.default,
        'FloatingMeltRate': basalforcings.default,
        'FrictionWaterPressure': friction.schoof,
        'FrontalForcingsRignotarma': frontalforcings.rignotarma,
        'FrontalForcingsSubglacialDischargearma': frontalforcings.rignotarma,
        'SMBarma': smb.arma,
        'SMBforcing': smb.default,
    }
    return structure if return_dict else list(structure.keys())

def _resolve_field(md, field=None, fieldname=None):
    """
    Retrieve a field either directly or via a dotted/indexed path.
    """

    if field is not None:
        return field, fieldname or "anonymous"

    if not fieldname:
        raise ValueError("Must specify either `field` or `fieldname`.")

    # Parse things like 'results.TransientSolution[0].Vel'
    attr_path = re.split(r'\[(.*?)\]', fieldname)[0]
    indexes = re.findall(r'\[(.*?)\]', fieldname)
    value = attrgetter(attr_path)(md)
    for idx in indexes:
        idx = idx.strip("'\" ")
        value = value[int(idx)] if idx.isdigit() else value[idx]
    return value, fieldname


def _check_size(md, field, fieldname, expected, message=None):
    """
    Check array size or 'universal' size semantics.
    """

    if isinstance(expected, str) and expected == "universal":
        v, e = md.mesh.numberofvertices, md.mesh.numberofelements
        if v == e or v + 1 == e or v == e + 1:
            md.check_message(message or f"{fieldname}: md.mesh.numberofvertices and md.mesh.numberofelements are within 1 of each other, cannot determine universal size")
        return

    if expected is not None:
        shape = np.shape(field)
        if len(shape) != len(expected) or any(
            not np.isnan(s) and s != shape[i] for i, s in enumerate(expected)
        ):
            md.check_message(message or f"{fieldname} has shape {shape}, expected {expected}")


def _check_values(md, field, fieldname, allowed, message=None):
    """
    Check categorical values.
    """

    if not np.all(np.isin(field, allowed)):
        md.check_message(message or f"{fieldname} values not in {allowed}")


def _check_bound(md, field, fieldname, op, bound, allow_nan, message=None):
    """
    Generic bound check: <, <=, >, >=.
    """

    fn = {
        ">": np.greater,
        ">=": np.greater_equal,
        "<": np.less,
        "<=": np.less_equal,
    }[op]

    field = np.asarray(field)
    
    # If allow_nan = True, only apply check on non-NaN values
    # (if all values are NaN, the check passes vacuously)
    if allow_nan:
        mask = ~np.isnan(field)
        allowed = np.all(fn(field[mask], bound))
    else:
    # allow_nan = False, apply check on all values
        allowed = np.all(fn(field, bound))
    
    if not allowed:
        md.check_message(message or f"{fieldname} fails condition {op} {bound}")

def _split_timeseries_field(md, field, kind):
    """
    Remove time row from timeseries field if present.
    """

    if field.ndim <= 1:
        return field, None

    nrow = field.shape[0]
    v, e = md.mesh.numberofvertices, md.mesh.numberofelements

    if kind == "timeseries":
        if nrow in (v + 1, e + 1):
            return field[:-1], field[-1]
        return field, None
    
    if kind == "singletimeseries":
        if nrow == 2:
            return field[0], field[1]
        return field, None

    if kind == "mappedtimeseries":
            return field[:-1], field[-1]
        
    return field, None


def _check_timeseries(md, field, fieldname, kind, message=None):
    """
    Check time series structure and sorting.
    """

    nrow = field.shape[0]
    v, e = md.mesh.numberofvertices, md.mesh.numberofelements

    def sorted_check(arr):
        if np.any(arr[:-1] > arr[1:]):
            md.check_message(message or f"{fieldname}: time not sorted")
        if np.any(arr[:-1] == arr[1:]):
            md.check_message(message or f"{fieldname}: duplicate timesteps")

    if kind == "timeseries":
        valid_rows = {v, e, v + 1, e + 1}
        if nrow not in valid_rows:
            md.check_message(message or f"{fieldname}: invalid timeseries row count")
        elif nrow in {v + 1, e + 1} and np.ndim(field) > 1:
            sorted_check(field[-1, :])

    elif kind == "singletimeseries":
        if nrow == 2:
            sorted_check(field[-1, :])
        elif nrow != 1:
            md.check_message(message or f"{fieldname}: must have 1 or 2 rows")

    elif kind == "mappedtimeseries":
        if np.ndim(field) > 1:
            sorted_check(field[-1, :])


def check_field(
    md,
    field = None,
    fieldname = None,
    allow_nan = True,
    allow_inf = True,
    scalar = False,
    size = None,
    numel = None,
    cell = False,
    string_list = False,
    allow_empty = False,
    values = None,
    gt = None,
    ge = None,
    lt = None,
    le = None,
    timeseries = False,
    singletimeseries = False,
    mappedtimeseries = False,
    file = False,
    message = None,
):

    """
    Validate a model field against a set of common checks.

    This function resolves a field either from the provided `field` argument
    or by resolving `fieldname` on `md` (supports dotted attribute paths and
    indexed access like "results.TransientSolution[0].Vel"). It then applies a
    series of optional validations (size, scalar, numel, emptiness, allowed
    values, numeric bounds, NaN/Inf presence, cell type, file existence and
    timeseries-specific structure checks). On failure, the model's
    md.check_message(...) method is invoked with an explanatory message.

    Parameters
    ----------
    md : :class:`object`
        ISSM model object. Must provide md.mesh.numberofvertices,
        md.mesh.numberofelements (for timeseries semantics) and a
        md.check_message(str) method to report validation failures.
    field : any, optional
        Direct field value to validate. If provided, ``fieldname`` is optional.
    fieldname : :class:`str`, optional
        Path used to resolve the field on ``md`` when ``field`` is ``None``.
        Supports dotted attribute access and index notation, e.g.
        "results.TransientSolution[0].Vel".
    allow_nan : :class:`bool`, optional
        If False, fail when the field contains NaNs. Default True.
    allow_inf : :class:`bool`, optional
        If False, fail when the field contains Infs. Default True.
    scalar : :class:`bool`, optional
        If True, require the field to be a scalar (or single element). Default False.
    size : :class:`tuple` or "universal", optional
        Expected shape of the field. Use np.nan for wildcard dimensions.
        The special string "universal" applies mesh-universal heuristics.
    numel : :class:`int` or sequence of int, optional
        Expected number of elements (np.size). Can be an int or list/tuple of
        acceptable ints.
    cell : :class:`bool`, optional
        If True, require the field to be a Python container (list/tuple/dict).
    string_list : :class:`bool`, optional
        If True, require the field to be a list of strings.
    allow_empty : :class:`bool`, optional
        If True, allow an empty field; otherwise an empty field triggers a check.
    values : ``sequence``, optional
        Allowed categorical values for the field (checked with np.isin).
    gt, ge, lt, le : :class:`scalar` or ``array-like``, optional
        Bound constraints (>, >=, <, <=). If provided, the corresponding
        comparison is applied elementwise.
    timeseries : :class:`bool`, optional
        Apply timeseries structural checks ("timeseries" kind).
    singletimeseries : :class:`bool`, optional
        Apply single-timeseries structural checks ("singletimeseries" kind).
    mappedtimeseries : :class:`bool`, optional
        Apply mapped-timeseries structural checks ("mappedtimeseries" kind).
    file : :class:`bool`, optional
        If True, treat the field as a filesystem path and require it exists.
    message : :class:`str`, optional
        Custom message to use on check failures instead of generated defaults.

    Returns
    -------
    md : object
        The same model object passed in (returned for convenience).

    Notes
    -----

    - If a resolved field is a Python scalar (:class:`bool`, :class:`int`, :class:`float`) it will be converted to a 1-element numpy array for uniform checking.
    - Checks do not raise exceptions; they report failures via md.check_message(msg).
    - The `size` parameter may contain np.nan entries that act as wildcards.
    """
    
    # Resolve field
    field, fieldname = _resolve_field(md, field, fieldname)

    if isinstance(field, (bool, int, float)):
        field = np.array([field])

    # Raw field checks
    # --------------------------------------------------------------
    # These checks are applied to the raw field as provided or resolved.

    # Empty
    if allow_empty and (field is None or len(np.atleast_1d(field)) == 0):
        md.check_message(message or f"{fieldname} is empty")

    # Size
    if size is not None:
        _check_size(md, field, fieldname, size, message)

    # Scalar
    if scalar:
        if not np.isscalar(field) and np.size(field) != 1:
            md.check_message(message or f"{fieldname} is not a scalar")

    # Numel
    if numel is not None:
        n = np.size(field)
        valid = [numel] if isinstance(numel, int) else numel
        if n not in valid:
            md.check_message(message or f"{fieldname} size {n}, expected {valid}")

    # NaN
    if not allow_nan and np.any(np.isnan(field)):
        md.check_message(message or f"{fieldname} contains NaN")

    # Inf
    if not allow_inf and np.any(np.isinf(field)):
        md.check_message(message or f"{fieldname} contains Inf")

    # Type
    if cell and not isinstance(field, (list, tuple, dict)):
        md.check_message(message or f"{fieldname} must be list/tuple/dict")

    # String_list
    if string_list and not isinstance(field, list):
        md.check_message(message or f"{fieldname} must be a list (string_list expected)")

    # Allowed values
    if values is not None:
        _check_values(md, field, fieldname, values, message)

    # File check
    if file and not os.path.exists(field):
        md.check_message(f"File in {fieldname} not found: {field}")

    # Timeseries checks
    if timeseries:
        _check_timeseries(md, field, fieldname, "timeseries", message)
    if singletimeseries:
        _check_timeseries(md, field, fieldname, "singletimeseries", message)
    if mappedtimeseries:
        _check_timeseries(md, field, fieldname, "mappedtimeseries", message)

    # Data-only field checks
    # --------------------------------------------------------------
    # When fields contain timeseries rows, the final row (time) is removed
    # before applying these checks.

    data_field = field
    kind = ("timeseries" if timeseries else
            "singletimeseries" if singletimeseries else
            "mappedtimeseries" if mappedtimeseries else None)
    if kind:
        data_field, _ = _split_timeseries_field(md, field, kind)
    
    # Bounds
    for op, bound in {">": gt, ">=": ge, "<": lt, "<=": le}.items():
        if bound is not None:
            _check_bound(md, data_field, fieldname, op, bound, allow_nan, message)

    return md

def cluster_queue_requirements(queue_dict, queue, np, time):
    """
    Validate cluster queue requirements.

    Validates that the requested queue exists, the time is positive and does not
    exceed the queue's maximum allowed time, and the number of processors is
    positive and does not exceed the queue's maximum allowed processors.

    Parameters
    ----------
    queue_dict : :class:`dict`
        Dictionary mapping queue names to (max_time, max_np) tuples, where
        max_time is the maximum allowed time for the queue and max_np is the
        maximum allowed number of processors.
    queue : :class:`str`
        Name of the requested queue.
    np : :class:`int`
        Number of processors requested. Must be positive and not exceed the
        queue's maximum.
    time : :class:`float` or :class:`int`
        Requested time in minutes (or appropriate time unit). Must be positive
        and not exceed the queue's maximum allowed time.

    Raises
    ------
    Exception

        - If the queue name is not found in queue_dict.
        - If time is not positive.
        - If requested time exceeds the maximum allowed time for the queue.
        - If number of processors is not positive.
        - If requested number of processors exceeds the maximum allowed for the queue.
    """

    try:
        rtime = queue_dict[queue][0]
    except KeyError:
        raise Exception(f'pyissm.classes.class_utils.cluster_queue_requirements: queue {queue} not recognized. Available queues are: {list(queue_dict.keys())}')
    
    if time <= 0:
        raise Exception('pyissm.classes.class_utils.cluster_queue_requirements: time must be positive')
    
    if time > rtime:
        raise Exception(f'pyissm.classes.class_utils.cluster_queue_requirements: requested time {time} exceeds maximum allowed time {rtime} for queue {queue}')
    
    if np <= 0:
        raise Exception('pyissm.classes.class_utils.cluster_queue_requirements: number of processors must be positive')

    max_np = queue_dict[queue][1]
    if np > max_np:
        raise Exception(f'pyissm.classes.class_utils.cluster_queue_requirements: requested number of processors {np} exceeds maximum allowed {max_np} for queue {queue}')