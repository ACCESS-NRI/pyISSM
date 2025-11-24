import numpy as np

CLASS_REGISTRY = {}

def register_class(cls):
    """
    Register a class in the global class registry for ISSM classes system.

    This decorator function automatically registers classes in the global CLASS_REGISTRY
    dictionary using a hierarchical naming scheme based on the module structure.
    It extracts the module path and creates a key for dynamic class instantiation.

    Parameters
    ----------
    cls : type
        The class to be registered in the registry.

    Returns
    -------
    type
        The same class that was passed in (decorator pattern).

    Notes
    -----
    The registration key is constructed as:
    - If 'classes' is in the module path: uses parts after 'classes' + class name
    - Otherwise: uses last 2 module parts + class name
    
    Example module path: pyissm.model.classes.smb -> key: "smb.default" (for class 'default')

    Examples
    --------
    @register_class
    class default(manage_state):
        pass
    # Registers as 'smb.default' if in pyissm.model.classes.smb module
    """
    parts = cls.__module__.split(".")
    if "classes" in parts:
        classes_index = parts.index("classes")
        key_parts = parts[classes_index + 1:]  # skip 'classes'
    else:
        key_parts = parts[-2:]  # fallback

    key = ".".join(key_parts + [cls.__name__])
    CLASS_REGISTRY[key] = cls
    return cls

def map_classtype(classtype: str) -> str:
    """
    Map legacy class types to modern hierarchical class names for backwards compatibility.

    This function provides a mapping from old ISSM class naming conventions to the
    new modular naming scheme. It ensures that legacy code and data files can still
    work with the updated class structure.

    Parameters
    ----------
    classtype : str
        The legacy class type string to be mapped.

    Returns
    -------
    str
        The modern class type string, or the original if no mapping exists.

    Notes
    -----
    The mapping covers major ISSM component categories:
    - basalforcings: Various basal forcing models (default, pico, linear, etc.)
    - calving: Ice calving models (default, crevassedepth, levermann, etc.)  
    - friction: Basal friction laws (default, coulomb, weertman, etc.)
    - smb: Surface mass balance models (default, arma, pdd, etc.)
    - hydrology: Subglacial hydrology models (glads, shreve, pism, etc.)
    - And many more...

    When a legacy class type is mapped, an informational message is printed.

    Examples
    --------
    >>> map_classtype('SMBforcing.SMBforcing')
    'ℹ️ Legacy classtype 'SMBforcing.SMBforcing' mapped to 'smb.default''
    'smb.default'
    
    >>> map_classtype('friction.weertman')  # Already modern
    'friction.weertman'
    """
    legacy_class_map = {
        'basalforcings.basalforcings': 'basalforcings.default',
        'basalforcingspico.basalforcingspico': 'basalforcings.pico',
        'linearbasalforcings.linearbasalforcings': 'basalforcings.linear',
        'linearbasalforcingsarma.linearbasalforcingsarma': 'basalforcings.lineararma',
        'mismipbasalforcings.mismipbasalforcings': 'basalforcings.mismip',
        'plumebasalforcings.plumebasalforcings': 'basalforcings.plume',
        'calving.calving': 'calving.default',
        'calvingcrevassedepth.calvingcrevassedepth': 'calving.crevassedepth',
        'calvingdev.calvingdev': 'calving.dev',
        'calvinglevermann.calvinglevermann': 'calving.levermann',
        'calvingparameterization.calvingparameterization': 'calving.parameterization',
        'calvingvonmises.calvingvonmises': 'calving.vonmises',
        'dsl.dsl': 'dsl.default',
        'dslmme.dslmme': 'dsl.mme',
        'love.love': 'love.default',
        'fourierlove': 'love.fourier',
        'friction.friction': 'friction.default',
        'frictioncoulomb.frictioncoulomb': 'friction.coulomb',
        'frictioncoulomb2.frictioncoulomb2': 'friction.coulomb2',
        'frictionhydro.frictionhydro': 'friction.hydro',
        'frictionjosh.frictionjosh': 'friction.josh',
        'frictionpism.frictionpism': 'friction.pism',
        'frictionregcoulomb.frictionregcoulomb': 'friction.regcoulomb',
        'frictionregcoulomb2.frictionregcoulomb2': 'friction.regcoulomb2',
        'frictionschoof.frictionschoof': 'friction.schoof',
        'frictionshakti.frictionshakti': 'friction.shakti',
        'frictionwaterlayer.frictionwaterlayer': 'friction.waterlayer',
        'frictionweertman.frictionweertman': 'friction.weertman',
        'frontalforcings.frontalforcings': 'frontalforcings.default',
        'frontalforcingsrignot.frontalforcingsrignot': 'frontalforcings.rignot',
        'frontalforcingsrignotarma.frontalforcingsrignotarma': 'frontalforcings.rignotarma',
        'generic.generic': 'cluster.generic',
        'hydrologyarmapw.hydrologyarmapw': 'hydrology.armapw',
        'hydrologydc.hydrologydc': 'hydrology.dc',
        'hydrologyglads.hydrologyglads': 'hydrology.glads',
        'hydrologypism.hydrologypism': 'hydrology.pism',
        'hydrologyshakti.hydrologyshakti': 'hydrology.shakti',
        'hydrologyshreve.hydrologyshreve': 'hydrology.shreve',
        'hydrologytws.hydrologytws': 'hydrology.tws',
        'inversion.inversion': 'inversion.default',
        'taoinversion.taoinversion': 'inversion.tao',
        'm1qn3inversion.m1qn3inversion': 'inversion.m1qn3',
        'matice.matice': 'materials.ice',
        'matdamageice.matdamageice': 'materials.damageice',
        'matenhancedice.matenhancedice': 'materials.enhancedice',
        'matestar.matestar': 'materials.estar',
        'mesh2d.mesh2d': 'mesh.mesh2d',
        'mesh3d.mesh3d': 'mesh.mesh3d',
        'mesh3dprisms.mesh3dprisms': 'mesh.mesh3dprisms',
        'mesh2dvertical.mesh2dvertical': 'mesh.mesh2dvertical',
        'mesh3dsurface.mesh3dsurface': 'mesh.mesh3dsurface',
        'SMBforcing.SMBforcing': 'smb.default',
        'SMBarma.SMBarma': 'smb.arma',
        'SMBcomponents.SMBcomponents': 'smb.components',
        'SMBd18opdd.SMBd18opdd': 'smb.d18opdd',
        'SMBgradients.SMBgradients': 'smb.gradients',
        'SMBgradientscomponents.SMBgradientscomponents': 'smb.gradientscomponents',
        'SMBgradientsela.SMBgradientsela': 'smb.gradientsela',
        'SMBhenning.SMBhenning': 'smb.henning',
        'SMBmeltcomponents.SMBmeltcomponents': 'smb.meltcomponents',
        'SMBpdd.SMBpdd': 'smb.pdd',
        'SMBpddSicopolis.SMBpddSicopolis': 'smb.pddSicopolis',
        'solidearth.solidearth': 'solidearth.earth', # TODO: Check for earth or europa?
        'solidearthsettings.solidearthsettings': 'solidearth.settings',
        'solidearthsolution.solidearthsolution': 'solidearth.solution',
        'timestepping.timestepping': 'timestepping.default',
        'timesteppingadaptive.timesteppingadaptive': 'timestepping.adaptive'
    }
    if classtype in legacy_class_map:
        print(f"ℹ️ Legacy classtype '{classtype}' mapped to '{legacy_class_map[classtype]}'")
    return legacy_class_map.get(classtype, classtype)


def create_instance(classtype: str):
    """
    Create an instance of a registered class by its type string.

    This function looks up a class in the CLASS_REGISTRY using the provided
    classtype string and instantiates it. It handles legacy class name mapping
    and provides error handling for unknown class types.

    Parameters
    ----------
    classtype : str
        The class type string identifying which class to instantiate.
        Can be either legacy or modern naming convention.

    Returns
    -------
    object or None
        An instance of the requested class, or None if the class type is unknown.

    Notes
    -----
    The function performs the following steps:
    1. Maps legacy class types to modern equivalents using map_classtype()
    2. Looks up the class in CLASS_REGISTRY
    3. Instantiates the class with no arguments
    4. Returns None (with warning) if class type is unknown

    This is the primary interface for dynamic class instantiation in the
    ISSM classes system, used for loading saved models and creating class instances
    from string identifiers.

    Examples
    --------
    >>> instance = create_instance('smb.default')
    >>> type(instance).__name__
    'default'
    
    >>> instance = create_instance('unknown.class')
    ⚠️ Unknown classtype unknown.class. Skipping...
    >>> instance is None
    True
    """
    classtype = map_classtype(classtype)
    if classtype not in CLASS_REGISTRY:
        print(f"⚠️ Unknown classtype {classtype}. Skipping...")
        return None
        # raise ValueError(f'create_instance: Unknown class type: {classtype}')
    return CLASS_REGISTRY[classtype]()

## Manage state for save/load and inheritance
class manage_state:
    """
    Base class providing state management and inheritance capabilities for ISSM classes.

    This class serves as the foundation for all ISSM classes, providing
    automatic field inheritance, state serialization for save/load operations,
    and field comparison utilities. It enables the "other" parameter pattern
    used throughout ISSM for class initialization and configuration inheritance.

    Methods
    -------
    __init__(self, other=None)
        Initialize the instance, optionally inheriting fields from another instance.
    _fields_equal(self, a, b)
        Compare two fields for equality, handling special cases like NaN and arrays.
    __getstate__(self)
        Get the current state for serialization (pickle support).
    __setstate__(self, state)
        Set the state from serialization (pickle support).

    Notes
    -----
    The inheritance mechanism works by:
    1. Comparing each attribute of the current instance with the "other" instance
    2. If the attribute exists in both and the values differ, the "other" value is used
    3. This allows default values to be overridden selectively

    Field equality comparison handles:
    - NaN values (considered equal to other NaN values)
    - NumPy arrays (using np.array_equal)
    - Regular Python objects (using == operator)

    Examples
    --------
    >>> class Example(manage_state):
    ...     def __init__(self, other=None):
    ...         self.value1 = 10
    ...         self.value2 = np.nan
    ...         super().__init__(other)
    
    >>> base = Example()
    >>> base.value1 = 20
    >>> inherited = Example(other=base)
    >>> inherited.value1
    20
    """

    ## Allow inheritance from existing model class
    def __init__(self, other=None):
        """
        Initialize the instance with optional field inheritance.

        Parameters
        ----------
        other : object, optional
            Another instance to inherit field values from. If provided, any
            fields in the current instance that differ from the default values
            and exist in 'other' will be replaced with the values from 'other'.

        Notes
        -----
        The inheritance process:
        1. Iterates through all attributes of the current instance
        2. Checks if the same attribute exists in the 'other' instance
        3. Compares field values using _fields_equal()
        4. If different, replaces the current value with the 'other' value

        This enables the common ISSM pattern where classes can inherit
        configurations from existing instances while maintaining their
        default values for unspecified fields.
        """

        # If other is provided...
        if other is not None:
            # Loop through all attributes of the current class...
            for attr in vars(self):
                # If the same attribute existing in the provided class, get the two fields
                if hasattr(other, attr):
                    field_other = getattr(other, attr)
                    field_self = getattr(self, attr)

                    # If the fields are different, replace the field_self with field_other
                    if not self._fields_equal(field_self, field_other):
                        setattr(self, attr, field_other)

    ## Check if fields are equal (used above)
    def _fields_equal(self, a, b):
        """
        Compare two field values for equality with special handling for NaN and arrays.

        Parameters
        ----------
        a : any
            First field value to compare.
        b : any  
            Second field value to compare.

        Returns
        -------
        bool
            True if the fields are considered equal, False otherwise.

        Notes
        -----
        Special cases handled:
        - NaN values: Two NaN values are considered equal
        - NumPy arrays: Compared using np.array_equal()
        - Other types: Compared using standard == operator

        This method is essential for the inheritance mechanism to work correctly
        with the variety of data types used in ISSM (scalars, arrays, NaN defaults).

        Examples
        --------
        >>> manager = manage_state()
        >>> manager._fields_equal(np.nan, np.nan)
        True
        >>> manager._fields_equal(np.array([1, 2]), np.array([1, 2]))
        True
        >>> manager._fields_equal(5, 5)
        True
        >>> manager._fields_equal(5, 10)
        False
        """

        # If both values are NaN...
        if isinstance(a, float) and np.isnan(a):
            return isinstance(b, float) and np.isnan(b)
        if isinstance(b, float) and np.isnan(b):
            return isinstance(a, float) and np.isnan(a)

        # If both values are equal arrays...
        if isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
            return np.array_equal(a, b)

        # If both values are equal scalars...
        return a == b

    ## Get the current state of self
    def __getstate__(self):
        """
        Get the current state of the instance for serialization.

        Returns
        -------
        dict
            Dictionary containing all instance attributes, suitable for pickling.

        Notes
        -----
        This method enables pickle support for ISSM classes, allowing models
        to be saved and loaded. It returns the instance's __dict__ which contains
        all the field values.

        Used by Python's pickle module and ISSM's save/load functionality.
        """
        return self.__dict__

    ## Set the current state of self
    def __setstate__(self, state):
        """
        Set the state of the instance from serialized data.

        Parameters
        ----------
        state : dict
            Dictionary containing the field values to restore.

        Notes
        -----
        This method enables pickle support for ISSM classes, allowing models
        to be loaded from saved files. It updates the instance's __dict__ with
        the provided state dictionary.

        Used by Python's pickle module and ISSM's save/load functionality.
        """
        self.__dict__.update(state)