from pyissm.model.classes import class_utils, class_registry
from pyissm.model import execute

@class_registry.register_class
class constants(class_registry.manage_state):
    """
    Physical constants class for ISSM.

    This class contains fundamental physical constants used in the ISSM framework.
    It provides standardized values for gravitational acceleration, Earth's rotation,
    time conversions, and other physical constants required for ice sheet modeling calculations.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    g : :class:`float`, default=9.81
        Gravitational acceleration [m/s^2].
    omega : :class:`float`, default=7.292e-5
        Angular velocity of Earth [rad/s].
    yts : :class:`float`, default=31536000.0
        Number of seconds in a year [s/yr] (365.0 * 24.0 * 3600.0).
    referencetemperature : :class:`float`, default=223.15
        Reference temperature used in the enthalpy model [K].
    gravitational_constant : :class:`float`, default=6.67259e-11
        Newtonian constant of gravitation [m^3/kg/s^2].

    Examples
    --------
    .. code-block:: python

        >>> md.constants = pyissm.model.classes.constants()
        >>> md.constants.g = 9.81
        >>> md.constants.yts = 365.25 * 24.0 * 3600.0
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.g = 9.81
        self.omega = 7.292 * 1e-5
        self.yts = 365.0 * 24.0 * 3600.0
        self.referencetemperature = 223.15
        self.gravitational_constant = 6.67259e-11

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   constants parameters:\n'
        s += '{}\n'.format(class_utils._field_display(self, 'g', 'gravitational acceleration [m/s^2]'))
        s += '{}\n'.format(class_utils._field_display(self, 'omega', 'angular velocity of Earth [rad/s]'))
        s += '{}\n'.format(class_utils._field_display(self, 'yts', 'number of seconds in a year [s/yr]'))
        s += '{}\n'.format(class_utils._field_display(self, 'referencetemperature', 'reference temperature used in the enthalpy model [K]'))
        s += '{}\n'.format(class_utils._field_display(self, 'gravitational_constant', 'Newtonian constant of gravitation [m^3/kg/s^2]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - constants Class'
        return s

    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [constants] parameters.

        Parameters
        ----------
        md : :class:`pyissm.model.Model`
            The model object to check.
        solution : :class:`str`
            The solution name to check.
        analyses : list of :class:`str`
            List of analyses to check consistency for.

        Returns 
        -------
        md : :class:`pyissm.model.Model`
            The model object with any consistency errors noted.
        """

        class_utils._check_field(md, fieldname = 'constants.g', ge = 0, scalar = True) # We allow 0 for validation tests
        class_utils._check_field(md, fieldname = 'constants.omega', ge = 0, scalar = True)
        class_utils._check_field(md, fieldname = 'constants.yts', ge = 0, scalar = True)
        class_utils._check_field(md, fieldname = 'constants.referencetemperature', scalar = True)
        class_utils._check_field(md, fieldname = 'constants.gravitational_constant', scalar = True)

        return md

    # Marshall method for saving the constants parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [constants] parameters to a binary file.

        Parameters
        ----------
        fid : :class:`file object`
            The file object to write the binary data to.
        prefix : :class:`str`
            Prefix string used for data identification in the binary file.
        md : :class:`pyissm.model.Model`, optional
            ISSM model object needed in some cases.
            
        Returns
        -------
        None
        """
        
        ## Write each field to the file (all fields are of the same type/format)
        fieldnames = ['g', 'omega', 'yts', 'referencetemperature', 'gravitational_constant']
        for fieldname in fieldnames:
            execute._write_model_field(fid, prefix, obj = self, fieldname = fieldname, format = 'Double')

