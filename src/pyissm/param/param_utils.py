import numpy as np

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
    md : object
        The model object containing the field to display.
    name : str
        The name of the field attribute to extract and display.
    comment : str
        Descriptive comment explaining the field's purpose and units.

    Returns
    -------
    str
        Formatted string representation of the field suitable for display.

    Examples
    --------
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
    offset : str
        String of spaces used for indentation in the output.
    name : str
        The name of the field being displayed.
    field : any
        The field value to be formatted for display.
    comment : str
        Descriptive comment explaining the field's purpose and units.

    Returns
    -------
    str
        Formatted string representation of the field.

    Notes
    -----
    Handles the following data types:
    - str: Shows string value (truncated if > 30 characters)
    - int, float: Shows numeric value
    - np.ndarray: Shows array shape
    - bool: Shows True/False
    - dict: Shows nested dictionary structure
    - list, tuple: Shows container contents (up to 5 elements)
    - None: Shows "None"
    - Other types: Shows "not displayed"

    Examples
    --------
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
    offset : str
        String of spaces used for base indentation in the output.
    name : str
        The name of the dictionary field being displayed.
    field : dict
        The dictionary to be formatted for display.
    comment : str
        Descriptive comment explaining the dictionary's purpose.

    Returns
    -------
    str
        Formatted multi-line string representation of the dictionary structure.

    Notes
    -----
    Empty dictionaries are displayed as 'N/A'. Non-empty dictionaries show
    '{dictionary}' on the first line followed by indented key-value pairs.

    Examples
    --------
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
    offset : str
        String of spaces used for indentation in the output.
    name : str
        The name of the list/tuple field being displayed.
    field : list or tuple
        The list or tuple to be formatted for display.
    comment : str
        Descriptive comment explaining the container's purpose.

    Returns
    -------
    str
        Formatted string representation of the list or tuple.

    Notes
    -----
    Small containers show individual elements: [1, 2, 3] or ('a', 'b', 'c').
    Large containers show size: [100x1] or (50x1).
    Only simple data types (str, bool, int, float) are shown individually.

    Examples
    --------
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
    offset : str
        String of spaces used for indentation.
    name : str
        The field name (truncated if > 23 characters).
    characterization : str
        The string representation of the field value.
    comment : str or list of str
        Descriptive comment(s) explaining the field. Can be a single string
        or a list of strings for multi-line comments.

    Returns
    -------
    str
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
    int
        The length of the longest field name in the object's __dict__.

    Examples
    --------
    >>> class Example:
    ...     def __init__(self):
    ...         self.a = 1
    ...         self.very_long_field_name = 2
    ...         self.b = 3
    >>> obj = Example()
    >>> getlongestfieldname(obj)
    20
    
    Notes
    -----
    Only considers attributes in the object's __dict__, not inherited
    attributes or properties.
    """

    maxlength = 0
    for key in self.__dict__.keys():
        length = len(key)
        if length > maxlength:
            maxlength = length

    return maxlength