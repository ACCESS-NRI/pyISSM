"""
Calving classes for ISSM.
"""

import numpy as np
import warnings
from pyissm.model.classes import class_utils, class_registry
from pyissm.model import execute, mesh

## ------------------------------------------------------
## calving.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    """
    Default calving parameters class for ISSM.

    This class contains the default parameters for calving in the ISSM framework.
    It defines the calving rate parameter.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    calvingrate : :class:`numpy.ndarray`, default=np.nan
        Calving rate at given location [m/a].

    Examples
    --------
    .. code-block:: python

        >>> md.calving = pyissm.model.classes.calving.default()
        >>> md.calving.calvingrate = np.zeros(md.mesh.numberofvertices, )
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.calvingrate = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving parameters:\n'

        s += '{}\n'.format(class_utils._field_display(self, 'calvingrate', 'calving rate at given location [m/a]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.default Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [calving.default] fields to 3D
        """
        self.calvingrate = mesh._project_3d(md, vector = self.calvingrate, type = 'node')
            
        return self

    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [calving.default] parameters.

        Parameters
        ----------
        md : :class:`pyissm.model.Model`
            The model object to check.
        solution : :class:`pyissm.model.solution`
            The solution object to check.
        analyses : list of :class:`str`
            List of analyses to check consistency for.

        Returns 
        -------
        md : :class:`pyissm.model.Model`
            The model object with any consistency errors noted.
        """

        #Early return if not transient with moving front
        if solution != 'TransientSolution' or not md.transient.ismovingfront:
            return md

        class_utils._check_field(md, fieldname = 'calving.calvingrate', ge = 0, timeseries = True, allow_nan = False, allow_inf = False)
            
        return md

    # Marshall method for saving the calving.default parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [calving.default] parameters to a binary file.

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

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute._write_model_field(fid, prefix, name = 'md.calving.law', data = 1, format = 'Integer')

        ## Write fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'calvingrate', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts, scale = 1. / md.constants.yts)

## ------------------------------------------------------
## calving.crevassedepth
## ------------------------------------------------------
@class_registry.register_class
class crevassedepth(class_registry.manage_state):
    """
    Crevasse depth calving parameters class for ISSM.

    This class contains the parameters for the crevasse depth calving model in the ISSM framework.
    It defines parameters related to crevasse opening stress, threshold, and water height.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    crevasse_opening_stress : :class:`int`, default=1
        0: stress only in the ice-flow direction, 1: maximum principal stress.
    crevasse_threshold : :class:`float`, default=1.0
        Ratio of full thickness to calve (e.g. 0.75 is for 75% of the total ice thickness).
    water_height : :class:`float`, default=0.0
        Water height in the crevasse [m].

    Examples
    --------
    .. code-block:: python
    
        >>> md.calving = pyissm.model.classes.calving.crevassedepth()
        >>> md.calving.crevasse_opening_stress = 1
        >>> md.calving.crevasse_threshold = 0.75
        >>> md.calving.water_height = 10.0
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.crevasse_opening_stress = 1
        self.crevasse_threshold = 1.
        self.water_height = 0.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving Pi parameters:\n'

        s += '{}\n'.format(class_utils._field_display(self, 'crevasse_opening_stress', '0: stress only in the ice-flow direction, 1: max principal'))
        s += '{}\n'.format(class_utils._field_display(self, 'crevasse_threshold', 'ratio of full thickness to calve (e.g. 0.75 is for 75% of the total ice thickness)'))
        s += '{}\n'.format(class_utils._field_display(self, 'water_height', 'water height in the crevasse [m]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.crevassedepth Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [calving.crevassedepth] fields to 3D
        """
        warnings.warn('pyissm.model.classes.calving.crevassedepth._extrude: 3D extrusion not implemented for calving.crevassedepth. Returning unchanged (2D) calving fields.')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [calving.crevassedepth] parameters.

        Parameters
        ----------
        md : :class:`pyissm.model.Model`
            The model object to check.
        solution : :class:`pyissm.model.solution`
            The solution object to check.
        analyses : list of :class:`str`
            List of analyses to check consistency for.

        Returns 
        -------
        md : :class:`pyissm.model.Model`
            The model object with any consistency errors noted.
        """

        # Early return if not transient with moving front
        if solution != "TransientSolution" or not md.transient.ismovingfront:
            return md

        class_utils._check_field(md, fieldname = "calving.crevasse_opening_stress", scalar = True, values = [0, 1])
        class_utils._check_field(md, fieldname = "calving.crevasse_threshold", scalar = True, gt = 0.0, le = 1.0)
        class_utils._check_field(md, fieldname = "calving.water_height", timeseries = True, ge = 0, allow_nan = False)

        return md

    # Marshall method for saving the calving.crevassedepth parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [calving.crevassedepth] parameters to a binary file.

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

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute._write_model_field(fid, prefix, name = 'md.calving.law', data = 6, format = 'Integer')

        ## Write fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'crevasse_opening_stress', format = 'Integer')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'crevasse_threshold', format = 'Double')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'water_height', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

## ------------------------------------------------------
## calving.dev
## ------------------------------------------------------
@class_registry.register_class
class dev(class_registry.manage_state):
    """
    Development calving parameters class for ISSM.

    This class contains the parameters for the development calving model in the ISSM framework.
    It defines stress thresholds for grounded and floating ice.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    stress_threshold_groundedice : :class:`float`, default=1e6
        Maximum stress threshold applied to grounded ice [Pa].
    stress_threshold_floatingice : :class:`float`, default=150e3
        Maximum stress threshold applied to floating ice [Pa].

    Examples
    --------
    .. code-block:: python
    
        >>> md.calving = pyissm.model.classes.calving.dev()
        >>> md.calving.stress_threshold_groundedice = 2e6
        >>> md.calving.stress_threshold_floatingice = 200e3
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.stress_threshold_groundedice = 1e6
        self.stress_threshold_floatingice = 150e3

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving Pi parameters:\n'

        s += '{}\n'.format(class_utils._field_display(self, 'stress_threshold_groundedice', 'sigma_max applied to grounded ice only [Pa]'))
        s += '{}\n'.format(class_utils._field_display(self, 'stress_threshold_floatingice', 'sigma_max applied to floating ice only [Pa]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.dev Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [calving.dev] fields to 3D
        """
        warnings.warn('pyissm.model.classes.calving.dev._extrude: 3D extrusion not implemented for calving.dev. Returning unchanged (2D) calving fields.')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [calving.dev] parameters.

        Parameters
        ----------
        md : :class:`pyissm.model.Model`
            The model object to check.
        solution : :class:`pyissm.model.solution`
            The solution object to check.
        analyses : list of :class:`str`
            List of analyses to check consistency for.

        Returns 
        -------
        md : :class:`pyissm.model.Model`
            The model object with any consistency errors noted.
        """

        # Early return if not transient with moving front
        if solution != "TransientSolution" or not md.transient.ismovingfront:
            return md

        class_utils._check_field(md, fieldname = "calving.stress_threshold_groundedice", gt = 0, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = "calving.stress_threshold_floatingice", gt = 0, allow_nan = False, allow_inf = False)

        return md

    # Marshall method for saving the calving.dev parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [calving.dev] parameters to a binary file.

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

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute._write_model_field(fid, prefix, name = 'md.calving.law', data = 2, format = 'Integer')

        ## Write fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'stress_threshold_groundedice', format = 'DoubleMat', mattype = 1)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'stress_threshold_floatingice', format = 'DoubleMat', mattype = 1)

## ------------------------------------------------------
## calving.levermann
## ------------------------------------------------------
@class_registry.register_class
class levermann(class_registry.manage_state):
    """
    Levermann calving parameters class for ISSM.

    This class contains the parameters for the Levermann calving model in the ISSM framework.
    It defines the proportionality coefficient used in the Levermann calving law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    coeff : :class:`float`, default=2e13
        Proportionality coefficient in the Levermann calving model.

    Examples
    --------
    .. code-block:: python
    
        >>> md.calving = pyissm.model.classes.calving.levermann()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.coeff = 2e13

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving Levermann parameters:\n'

        s += '{}\n'.format(class_utils._field_display(self, 'coeff', 'proportionality coefficient in Levermann model'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.levermann Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [calving.levermann] fields to 3D
        """
        self.coeff = mesh._project_3d(md, vector = self.coeff, type = 'node')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [calving.levermann] parameters.

        Parameters
        ----------
        md : :class:`pyissm.model.Model`
            The model object to check.
        solution : :class:`pyissm.model.solution`
            The solution object to check.
        analyses : list of :class:`str`
            List of analyses to check consistency for.

        Returns 
        -------
        md : :class:`pyissm.model.Model`
            The model object with any consistency errors noted.
        """

        # Early return if not transient with moving front
        if solution != "TransientSolution" or not md.transient.ismovingfront:
            return md

        class_utils._check_field(md, fieldname = "calving.coeff", size = (md.mesh.numberofvertices, ), gt = 0)

        return md

    # Marshall method for saving the calving.levermann parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [calving.levermann] parameters to a binary file.

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

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute._write_model_field(fid, prefix, name = 'md.calving.law', data = 3, format = 'Integer')

        ## Write fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'coeff', format = 'DoubleMat', mattype = 1)

## ------------------------------------------------------
## calving.minthickness
## ------------------------------------------------------
@class_registry.register_class
class minthickness(class_registry.manage_state):
    """
    Minimum thickness calving parameters class for ISSM.

    This class contains the parameters for the minimum thickness calving model in the ISSM framework.
    It defines the minimum ice thickness below which no ice is allowed.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    min_thickness : class:`float`, default=100.
        Minimum thickness below which no ice is allowed [m].

    Examples
    --------
    .. code-block:: python

        >>> md.calving = pyissm.model.classes.calving.minthickness()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.min_thickness = 100.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving Minimum thickness:\n'

        s += '{}\n'.format(class_utils._field_display(self, 'min_thickness', 'minimum thickness below which no ice is allowed'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.minthickness Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [calving.minthickness] fields to 3D
        """
        warnings.warn('pyissm.model.classes.calving.minthickness._extrude: 3D extrusion not implemented for calving.minthickness. Returning unchanged (2D) calving fields.')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [calving.minthickness] parameters.

        Parameters
        ----------
        md : :class:`pyissm.model.Model`
            The model object to check.
        solution : :class:`pyissm.model.solution`
            The solution object to check.
        analyses : list of :class:`str`
            List of analyses to check consistency for.

        Returns 
        -------
        md : :class:`pyissm.model.Model`
            The model object with any consistency errors noted.
        """

        # Early return if not transient with moving front
        if solution != "TransientSolution" or not md.transient.ismovingfront:
            return md

        class_utils._check_field(md, fieldname = "calving.min_thickness", gt = 0, scalar = True, allow_nan = False, allow_inf = False)

        return md
    
    # Marshall method for saving the calving.minthickness parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [calving.minthickness] parameters to a binary file.

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

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute._write_model_field(fid, prefix, name = 'md.calving.law', data = 4, format = 'Integer')

        ## Write fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'min_thickness', format = 'Double')

## ------------------------------------------------------
## calving.parameterization
## ------------------------------------------------------
@class_registry.register_class
class parameterization(class_registry.manage_state):
    """
    Parameterization calving parameters class for ISSM.

    This class contains the parameters for the parameterization calving model in the ISSM framework.
    It defines parameters controlling the calving rate as a function of ice thickness, velocity, and other factors.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    min_thickness : :class:`float`, default=0.
        Minimum thickness below which no ice is allowed [m].
    use_param : :class:`float`, default=0.
        Parameterization mode selector:
            -1: use frontal ablation rate,
             0: linear function,
             1: tanh function,
             2: tanh(thickness),
             3: tanh(normalized velocity),
             4: tanh(truncated velocity),
             5: linear(truncated velocity).
    theta : :class:`float`, default=0.
        Amplifier parameter for the calving law.
    alpha : :class:`float`, default=0.
        Slope parameter for the calving law.
    xoffset : :class:`float`, default=0.
        Offset in x-axis for the calving law.
    yoffset : :class:`float`, default=0.
        Offset in y-axis for the calving law.
    vel_upperbound : :class:`float`, default=6000.
        Upper bound of ice velocity to reduce the calving rate [m/a].
    vel_threshold : :class:`float`, default=0.
        Threshold of ice velocity to reduce the calving rate [m/a].
    vel_lowerbound : :class:`float`, default=0.
        Lower bound of ice velocity to reduce the calving rate [m/a].

    Examples
    --------
    .. code-block:: python
    
        >>> md.calving = pyissm.model.classes.calving.parameterization()
        >>> md.calving.min_thickness = 50.
        >>> md.calving.use_param = 1
        >>> md.calving.theta = 2.0
        >>> md.calving.alpha = 0.01
        >>> md.calving.xoffset = 10.0
        >>> md.calving.yoffset = 0.5
        >>> md.calving.vel_upperbound = 5000.
        >>> md.calving.vel_threshold = 1000.
        >>> md.calving.vel_lowerbound = 100.
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.min_thickness = 0.
        self.use_param = 0.
        self.theta = 0.
        self.alpha = 0.
        self.xoffset = 0.
        self.yoffset = 0.
        self.vel_upperbound = 6000.
        self.vel_threshold = 0.
        self.vel_lowerbound = 0.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving test parameters:\n'

        s += '{}\n'.format(class_utils._field_display(self, 'min_thickness', 'minimum thickness below which no ice is allowed [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'use_param', '-1 - just use frontal ablation rate, 0 - f(x) = y_{o} + \alpha (x+x_{o}), 1 - f(x)=y_{o}-\frac{\theta}{2}\tanh(\alpha(x+x_{o})), 2 - tanh(thickness), 3 - tanh(normalized vel), 4 - tanh(truncated vel), 5 - linear(truncated vel)'))
        s += '{}\n'.format(class_utils._field_display(self, 'theta', 'the amplifier'))
        s += '{}\n'.format(class_utils._field_display(self, 'alpha', 'the slope'))
        s += '{}\n'.format(class_utils._field_display(self, 'xoffset', 'offset in x-axis'))
        s += '{}\n'.format(class_utils._field_display(self, 'yoffset', 'offset in y-axis'))
        s += '{}\n'.format(class_utils._field_display(self, 'vel_lowerbound', 'lowerbound of ice velocity to reduce the calving rate [m/a]'))
        s += '{}\n'.format(class_utils._field_display(self, 'vel_threshold', 'threshold of ice velocity to reduce the calving rate [m/a]'))
        s += '{}\n'.format(class_utils._field_display(self, 'vel_upperbound', 'upperbound of ice velocity to reduce the calving rate [m/a]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.parameterization Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [calving.parameterization] fields to 3D
        """
        warnings.warn('pyissm.model.classes.calving.parameterization._extrude: 3D extrusion not implemented for calving.parameterization. Returning unchanged (2D) calving fields.')
            
        return self

    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [calving.parameterization] parameters.

        Parameters
        ----------
        md : :class:`pyissm.model.Model`
            The model object to check.
        solution : :class:`pyissm.model.solution`
            The solution object to check.
        analyses : list of :class:`str`
            List of analyses to check consistency for.

        Returns 
        -------
        md : :class:`pyissm.model.Model`
            The model object with any consistency errors noted.
        """

        # Early return if not transient with moving front
        if solution != "TransientSolution" or not md.transient.ismovingfront:
            return md

        class_utils._check_field(md, fieldname = "calving.min_thickness", scalar = True, ge = 0, allow_nan = False)
        class_utils._check_field(md, fieldname = "calving.use_param", scalar = True, values = [-1, 0, 1, 2, 3, 4, 5])
        class_utils._check_field(md, fieldname = "calving.theta", scalar = True, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = "calving.alpha", scalar = True, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = "calving.xoffset", scalar = True, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = "calving.yoffset", scalar = True, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = "calving.vel_lowerbound", scalar = True, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = "calving.vel_threshold", scalar = True, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = "calving.vel_upperbound", scalar = True, allow_nan = False, allow_inf = False)

        return md
    
    # Marshall method for saving the calving.parameterization parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [calving.parametrization] parameters to a binary file.

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

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute._write_model_field(fid, prefix, name = 'md.calving.law', data = 9, format = 'Integer')

        ## Write Double fields
        fieldnames = ['min_thickness', 'use_param', 'theta', 'alpha', 'xoffset', 'yoffset']
        for field in fieldnames:
            execute._write_model_field(fid, prefix, obj = self, fieldname = field, format = 'Double')

        execute._write_model_field(fid, prefix, obj = self, fieldname = 'vel_lowerbound', format = 'Double', scale = 1. / md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'vel_threshold', format = 'Double', scale = 1. / md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'vel_upperbound', format = 'Double', scale = 1. / md.constants.yts)

## ------------------------------------------------------
## calving.vonmises
## ------------------------------------------------------
@class_registry.register_class
class vonmises(class_registry.manage_state):
    """
    Von Mises calving parameters class for ISSM.

    This class contains the parameters for the Von Mises calving model in the ISSM framework.
    It defines stress thresholds for grounded and floating ice, as well as the minimum ice thickness below which no ice is allowed.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    stress_threshold_groundedice : :class:`float`, default=0
        Maximum Von Mises stress threshold applied to grounded ice [Pa].
    stress_threshold_floatingice : :class:`float`, default=0
        Maximum Von Mises stress threshold applied to floating ice [Pa].
    min_thickness : :class:`float`, default=0.
        Minimum thickness below which no ice is allowed [m].
        
    Examples
    --------
    .. code-block:: python
    
        >>> md.calving = pyissm.model.classes.calving.vonmises()
        >>> md.calving.stress_threshold_groundedice = 1e6
        >>> md.calving.stress_threshold_floatingice = 150e3
        >>> md.calving.min_thickness = 50.
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.stress_threshold_groundedice = 1e6
        self.stress_threshold_floatingice = 150e3
        self.min_thickness = 0.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Calving VonMises parameters:\n'

        s += '{}\n'.format(class_utils._field_display(self, 'stress_threshold_groundedice', 'sigma_max applied to grounded ice only [Pa]'))
        s += '{}\n'.format(class_utils._field_display(self, 'stress_threshold_floatingice', 'sigma_max applied to floating ice only [Pa]'))
        s += '{}\n'.format(class_utils._field_display(self, 'min_thickness', 'minimum thickness below which no ice is allowed [m]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - calving.vonmises Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [calving.vonmises] fields to 3D
        """
        warnings.warn('pyissm.model.classes.calving.vonmises._extrude: 3D extrusion not implemented for calving.vonmises. Returning unchanged (2D) calving fields.')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [calving.vonmises] parameters.

        Parameters
        ----------
        md : :class:`pyissm.model.Model`
            The model object to check.
        solution : :class:`pyissm.model.solution`
            The solution object to check.
        analyses : list of :class:`str`
            List of analyses to check consistency for.

        Returns 
        -------
        md : :class:`pyissm.model.Model`
            The model object with any consistency errors noted.
        """

        # Early return if not transient with moving front
        if solution != "TransientSolution" or not md.transient.ismovingfront:
            return md

        class_utils._check_field(md, fieldname = "calving.stress_threshold_groundedice", gt = 0, size = 'universal', allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = "calving.stress_threshold_floatingice", gt = 0, size = 'universal', allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = "calving.min_thickness", scalar = True, ge = 0, allow_nan = False, allow_inf = False)

        return md
    
    # Marshall method for saving the calving.vonmises parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [calving.vonmises] parameters to a binary file.

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

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute._write_model_field(fid, prefix, name = 'md.calving.law', data = 2, format = 'Integer')

        ## Write fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'stress_threshold_groundedice', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'stress_threshold_floatingice', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'min_thickness', format = 'Double')