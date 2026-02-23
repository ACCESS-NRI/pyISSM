"""
Friction classes for ISSM.
"""

import numpy as np
from pyissm.model.classes import class_utils, class_registry
from pyissm.model import execute, mesh

## ------------------------------------------------------
## friction.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    """
    Default friction class for ISSM.

    This class contains the parameters for the default (Budd) friction law in the ISSM framework.
    It defines the main friction-related parameters specific to the Budd law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    coefficient : :class:`nump.ndarray`, default=np.nan
        Friction coefficient [SI].
    p : :class:`nump.ndarray`, default=np.nan
        p exponent.
    q : :class:`nump.ndarray`, default=np.nan
        q exponent.
    coupling : :class:`int`, default=0
        Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: used coupled model (not implemented yet).
    linearize : :class:`int`, default=0
        0: not linearized, 1: interpolated linearly, 2: constant per element (default is 0).
    effective_pressure : :class:`nump.ndarray`, default=np.nan
        Effective Pressure for the forcing if not coupled [Pa].
    effective_pressure_limit : :class:`nump.ndarray`, default=0
        Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0).

    Examples
    --------
    .. code-block:: python

        >>> md.friction = pyissm.model.classes.friction.default()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.coefficient = np.nan
        self.p = np.nan
        self.q = np.nan
        self.coupling = 0
        self.linearize = 0
        self.effective_pressure = np.nan
        self.effective_pressure_limit = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Basal shear stress parameters: Sigma_b = coefficient^2 * Neff ^r * |u_b|^(s - 1) * u_b,\n'
        s += '(effective stress Neff = rho_ice * g * thickness + rho_water * g * base, r = q / p and s = 1 / p)\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'coefficient', 'friction coefficient [SI]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'p', 'p exponent'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'q', 'q exponent'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'coupling', 'Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: used coupled model (not implemented yet)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'linearize', '0: not linearized, 1: interpolated linearly, 2: constant per element (default is 0)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.default Class'
        return s
    
    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude friction.default fields to 3D
        """
        self.coefficient = mesh.project_3d(md, vector = self.coefficient, type = 'node', layer = 1)
        self.p = mesh.project_3d(md, vector = self.p, type = 'element')
        self.q = mesh.project_3d(md, vector = self.q, type = 'element')
        if self.coupling in[3, 4]:
            self.effective_pressure = mesh.project3d(md, vector = self.effective_pressure, type = 'node', layer = 1)

        return self

    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [friction.default] parameters.

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

        # Early return if necessary analyses or solutions not specified
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        if solution == 'TransientSolution' and not md.transient.isstressbalance and not md.transient.isthermal:
            return md
        
        class_utils.check_field(md, fieldname = "friction.coefficient", timeseries = True, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "friction.p", size = (md.mesh.numberofelements, ), allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "friction.q", size = (md.mesh.numberofelements, ), allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "friction.coupling", scalar = True, values = [0, 1, 2, 3, 4])
        class_utils.check_field(md, fieldname = "friction.linearize", scalar = True, values = [0, 1, 2])
        class_utils.check_field(md, fieldname = "friction.effective_pressure_limit", scalar = True, ge = 0)

        if self.coupling == 3:
            class_utils.check_field(md, fieldname = "friction.effective_pressure", timeseries = True, allow_nan = False, allow_inf = False)

        return md

    # Marshall method for saving the friction.default parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.default] parameters to a binary file.

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

        ## Write headers to file
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, data = 1, name = 'md.friction.law', format = 'Integer')

        ## Write coefficient field
        if isinstance(self.coefficient, np.ndarray) and ((self.coefficient.shape[0] == md.mesh.numberofvertices) or (self.coefficient.shape[0] == md.mesh.numberofvertices + 1)):
            execute.WriteData(fid, prefix, obj = self, fieldname = 'coefficient', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        elif isinstance(self.coefficient, np.ndarray) and ((self.coefficient.shape[0] == md.mesh.numberofelements) or (self.coefficient.shape[0] == md.mesh.numberofelements + 1)):
            execute.WriteData(fid, prefix, obj = self, fieldname = 'coefficient', format = 'DoubleMat', mattype = 2, timeserieslength = md.mesh.numberofelements + 1, yts = md.constants.yts)
        else:
            raise RuntimeError('friction coefficient time series should be a vertex or element time series')
        
        ## Write other fields with specific formats
        execute.WriteData(fid, prefix, obj = self, fieldname = 'p', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'q', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'coupling', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'linearize', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'effective_pressure_limit', format = 'Double')
        
        ## Write conditional effective pressure
        if (self.coupling == 3) or (self.coupling == 4):
            execute.WriteData(fid, prefix, obj = self, fieldname = 'effective_pressure', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        elif self.coupling > 4:
            raise ValueError(f'md.friction.coupling = {self.coupling} is not implemented yet')
        

## ------------------------------------------------------
## friction.coulomb
## ------------------------------------------------------
@class_registry.register_class
class coulomb(class_registry.manage_state):
    """
    Coulomb friction class for ISSM.

    This class contains the parameters for the Coulomb friction law in the ISSM framework.
    It defines the main friction-related parameters specific to the Coulomb law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    coefficient : :class:`numpy.ndarray`, default=np.nan
        Power law (Weertman) friction coefficient [SI].
    coefficientcoulomb : :class:`numpy.ndarray`, default=np.nan
        Coulomb friction coefficient [SI].
    p : :class:`numpy.ndarray`, default=np.nan
        p exponent.
    q : :class:`numpy.ndarray`, default=np.nan
        q exponent.
    coupling : :class:`int`, default=0
        Coupling flag: 0 for default, 1 for forcing (provide md.friction.effective_pressure), 2 for coupled (not implemented yet).
    effective_pressure : :class:`numpy.ndarray`, default=np.nan
        Effective Pressure for the forcing if not coupled [Pa].
    effective_pressure_limit : :class:`numpy.ndarray`, default=0
        Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0).

    Examples
    --------
    .. code-block:: python
        
        >>> md.friction = pyissm.model.classes.friction.coulomb()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.coefficient = np.nan
        self.coefficientcoulomb = np.nan
        self.p = np.nan
        self.q = np.nan
        self.coupling = 0
        self.effective_pressure = np.nan
        self.effective_pressure_limit = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Basal shear stress parameters: Sigma_b = min(coefficient^2 * Neff ^r * |u_b|^(s - 1) * u_b,\n'
        s += 'coefficientcoulomb^2 * Neff), (effective stress Neff = rho_ice * g * thickness + rho_water * g * bed, r = q / p and s = 1 / p).\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'coefficient', 'power law (Weertman) friction coefficient [SI]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'coefficientcoulomb', 'Coulomb friction coefficient [SI]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'p', 'p exponent'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'q', 'q exponent'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'coupling', 'Coupling flag: 0 for default, 1 for forcing(provide md.friction.effective_pressure)  and 2 for coupled(not implemented yet)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.coulomb Class'
        return s
    
    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude friction.coulomb fields to 3D
        """
        self.coefficient = mesh.project_3d(md, vector = self.coefficient, type = 'node',  layer = 1)
        self.coefficientcoulomb = mesh.project_3d(md, vector = self.coefficientcoulomb, type = 'node', layer = 1)
        self.p = mesh.project_3d(md, vector = self.p, type = 'element')
        self.q = mesh.project_3d(md, vector = self.q, type = 'element')
        if self.coupling == 1:
            self.effective_pressure = mesh.project_3d(md, vector = self.effective_pressure, type = 'node', layer = 1)
        elif self.coupling >= 2:
            raise ValueError('pyissm.model.classes.friction.coulomb.extrude: md.friction.coupling >= 2 not implemented yet.')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [friction.coulomb] parameters.

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

        # Early return if necessary analyses or solutions not specified
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        
        class_utils.check_field(md, fieldname = "friction.coefficient", timeseries = True, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "friction.coefficientcoulomb", timeseries = True, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "friction.p", size = (md.mesh.numberofelements, ), allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "friction.q", size = (md.mesh.numberofelements, ), allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "friction.coupling", scalar = True, values = [0, 1, 2])
        class_utils.check_field(md, fieldname = "friction.effective_pressure_limit", scalar = True, ge = 0)

        if self.coupling == 1:
            class_utils.check_field(md, fieldname = "friction.effective_pressure", timeseries = True, allow_nan = False, allow_inf = False)
        elif self.coupling == 2:
            raise ValueError('pyissm.model.classes.friction.coulomb.check_consistency: md.friction.coupling = 2 (coupled) is not implemented yet')
        elif self.coupling > 2:
            raise ValueError(f'pyissm.model.classes.friction.coulomb.check_consistency: md.friction.coupling = {self.coupling} is not implemented yet')

        return md

    # Marshall method for saving the friction.coulomb parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.coulomb] parameters to a binary file.

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

        ## Write headers to file
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.friction.law', data = 7, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'coefficient', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'coefficientcoulomb', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'p', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'q', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'coupling', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'effective_pressure_limit', format = 'Double')

        ## Write conditional effective pressure
        if self.coupling == 1:
            execute.WriteData(fid, prefix, obj = self, fieldname = 'effective_pressure', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        elif self.coupling > 1:
            raise ValueError(f'md.friction.coupling = {self.coupling} is not implemented yet')

## ------------------------------------------------------
## friction.coulomb2
## ------------------------------------------------------
@class_registry.register_class
class coulomb2(class_registry.manage_state):
    """
    Coulomb2 friction class for ISSM.

    This class contains the parameters for the Coulomb2 friction law in the ISSM framework.
    It defines the main friction-related parameters specific to the Coulomb2 law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    coefficient : :class:`numpy.ndarray`, default=np.nan
        Power law (Weertman) friction coefficient [SI].
    coefficientcoulomb : :class:`numpy.ndarray`, default=np.nan
        Coulomb friction coefficient [SI].
    p : :class:`numpy.ndarray`, default=np.nan
        p exponent.
    q : :class:`numpy.ndarray`, default=np.nan
        q exponent.
    coupling : :class:`int`, default=0
        Coupling flag: 0 for default, 1 for forcing (provide md.friction.effective_pressure), 2 for coupled (not implemented yet).
    effective_pressure : :class:`numpy.ndarray`, default=np.nan
        Effective Pressure for the forcing if not coupled [Pa].
    effective_pressure_limit : :class:`numpy.ndarray`, default=0
        Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0).

    Examples
    --------
    .. code-block:: python
    
        >>> md.friction = pyissm.model.classes.friction.coulomb2()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.coefficient = np.nan
        self.coefficientcoulomb = np.nan
        self.p = np.nan
        self.q = np.nan
        self.coupling = 0
        self.effective_pressure = np.nan
        self.effective_pressure_limit = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Basal shear stress parameters: Sigma_b = min(coefficient^2 * Neff ^r * |u_b|^(s - 1) * u_b,\n'
        s += 'coefficientcoulomb^2 * Neff), (effective stress Neff = rho_ice * g * thickness + rho_water * g * bed, r = q / p and s = 1 / p).\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'coefficient', 'power law (Weertman) friction coefficient [SI]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'coefficientcoulomb', 'Coulomb friction coefficient [SI]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'p', 'p exponent'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'q', 'q exponent'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'coupling', 'Coupling flag: 0 for default, 1 for forcing(provide md.friction.effective_pressure)  and 2 for coupled(not implemented yet)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.coulomb2 Class'
        return s

    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude friction.coulomb2 fields to 3D
        """
        self.coefficient = mesh.project_3d(md, vector = self.coefficient, type = 'node',  layer = 1)
        self.coefficientcoulomb = mesh.project_3d(md, vector = self.coefficientcoulomb, type = 'node', layer = 1)
        self.p = mesh.project_3d(md, vector = self.p, type = 'element')
        self.q = mesh.project_3d(md, vector = self.q, type = 'element')
        if self.coupling == 1:
            self.effective_pressure = mesh.project_3d(md, vector = self.effective_pressure, type = 'node', layer = 1)
        elif self.coupling >= 2:
            raise ValueError('pyissm.model.classes.friction.coulomb2.extrude: md.friction.coupling >= 2 not implemented yet.')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [friction.coulomb2] parameters.

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

        # Early return if necessary analyses or solutions not specified
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        
        class_utils.check_field(md, fieldname = "friction.coefficient", timeseries = True, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "friction.coefficientcoulomb", timeseries = True, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "friction.p", size = (md.mesh.numberofelements, ), allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "friction.q", size = (md.mesh.numberofelements, ), allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "friction.coupling", scalar = True, values = [0, 1, 2])
        class_utils.check_field(md, fieldname = "friction.effective_pressure_limit", scalar = True, ge = 0)

        if self.coupling == 1:
            class_utils.check_field(md, fieldname = "friction.effective_pressure", timeseries = True, allow_nan = False, allow_inf = False)
        elif self.coupling == 2:
            raise ValueError('pyissm.model.classes.friction.coulomb.check_consistency: md.friction.coupling = 2 (coupled) is not implemented yet')
        elif self.coupling > 2:
            raise ValueError(f'pyissm.model.classes.friction.coulomb.check_consistency: md.friction.coupling = {self.coupling} is not implemented yet')

        return md

    # Marshall method for saving the friction.coulomb2 parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.coulomb2] parameters to a binary file.

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

        ## Write headers to file
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.friction.law', data = 7, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'coefficient', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'coefficientcoulomb', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'p', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'q', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'coupling', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'effective_pressure_limit', format = 'Double')

        ## Write conditional effective pressure
        if self.coupling == 1:
            execute.WriteData(fid, prefix, obj = self, fieldname = 'effective_pressure', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        elif self.coupling > 1:
            raise ValueError(f'md.friction.coupling = {self.coupling} is not implemented yet')

## ------------------------------------------------------
## friction.hydro
## ------------------------------------------------------
@class_registry.register_class
class hydro(class_registry.manage_state):
    """
    Hydro friction class for ISSM.

    This class contains the parameters for the hydro (Gagliardini 2007) friction law in the ISSM framework.
    It defines the main friction-related parameters specific to the hydro law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    coupling : :class:`int`, default=0
        Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: use coupled model (not implemented yet).
    q : :class:`numpy.ndarray`, default=np.nan
        Friction law exponent q >= 1.
    C : :class:`numpy.ndarray`, default=np.nan
        Friction law max value (Iken bound).
    As : :class:`numpy.ndarray`, default=np.nan
        Sliding parameter without cavitation [m Pa^-n s^-1].
    effective_pressure : :class:`numpy.ndarray`, default=np.nan
        Effective Pressure for the forcing if not coupled [Pa].
    effective_pressure_limit : :class:`numpy.ndarray`, default=0
        Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0).

    Examples
    --------
    .. code-block:: python

        >>> md.friction = pyissm.model.classes.friction.hydro()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.coupling = 0
        self.q = np.nan
        self.C = np.nan
        self.As = np.nan
        self.effective_pressure = np.nan
        self.effective_pressure_limit = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Effective Pressure based friction law described in Gagliardini 2007\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'coupling', 'Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: used coupled model (not implemented yet)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'q', 'friction law exponent q >= 1'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'C', 'friction law max value (Iken bound)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'As', 'Sliding Parameter without cavitation [m Pa^ - n s^ - 1]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.hydro Class'
        return s
    
    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude friction.hydro fields to 3D
        """
        self.q = mesh.project_3d(md, vector = self.q, type = 'element')
        self.C = mesh.project_3d(md, vector = self.C, type = 'element')
        self.As = mesh.project_3d(md, vector = self.As, type = 'element')
        self.q = mesh.project_3d(md, vector = self.q, type = 'element')
        if self.coupling in [3, 4]:
            self.effective_pressure = mesh.project_3d(md, vector = self.effective_pressure, type = 'node', layer = 1)
        elif self.coupling > 2:
            raise ValueError('pyissm.model.classes.friction.hydro.extrude: md.friction.coupling > 4 not implemented yet.')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [friction.hydro] parameters.

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

        # Early return if necessary analyses or solutions not specified
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        
        class_utils.check_field(md, fieldname = "friction.coupling", scalar = True, values = [0, 1, 2, 3, 4])
        class_utils.check_field(md, fieldname = "friction.q", size = (md.mesh.numberofelements, ), allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "friction.C", size = (md.mesh.numberofelements, ), allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "friction.As", size = (md.mesh.numberofelements, ), allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "friction.effective_pressure_limit", scalar = True, ge = 0)       

        if self.coupling == 3:
            class_utils.check_field(md, fieldname = "friction.effective_pressure", timeseries = True, allow_nan = False, allow_inf = False)
        elif self.coupling > 4:
            raise ValueError(f'pyissm.model.classes.friction.hydro.check_consistency: md.friction.coupling = {self.coupling} is not implemented yet. Use md.friction.coupling <= 4.')

        return md

    # Marshall method for saving the friction.hydro parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.hydro] parameters to a binary file.

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

        ## Write headers to file
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.friction.law', data = 3, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'q', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'C', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'As', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'coupling', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'effective_pressure_limit', format = 'Double')

        ## Write conditional effective pressure
        if self.coupling in [3, 4]:
            execute.WriteData(fid, prefix, obj = self, fieldname = 'effective_pressure', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        elif self.coupling > 4:
            raise ValueError(f'md.friction.coupling = {self.coupling} is not implemented yet')

## ------------------------------------------------------
## friction.josh
## ------------------------------------------------------
@class_registry.register_class
class josh(class_registry.manage_state):
    """
    Josh friction class for ISSM.

    This class contains the parameters for the Josh friction law in the ISSM framework.
    It defines the main friction-related parameters specific to the Josh law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    coefficient : :class:`numpy.ndarray`, default=np.nan
        Friction coefficient [SI].
    pressure_adjusted_temperature : :class:`numpy.ndarray`, default=np.nan
        Friction pressure_adjusted_temperature (T - Tpmp) [K].
    gamma : :class:`numpy.ndarray`, default=1.
        (T - Tpmp)/gamma [K].
    effective_pressure_limit : :class:`numpy.ndarray`, default=0
        Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0).

    Examples
    --------
    .. code-block:: python

        >>> md.friction = pyissm.model.classes.friction.josh()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.coefficient = np.nan
        self.pressure_adjusted_temperature = np.nan
        self.gamma = 1.
        self.effective_pressure_limit = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Basal shear stress parameters: Sigma_b = coefficient^2 * Neff ^r * |u_b|^(s - 1) * u_b,\n'
        s += '(effective stress Neff = rho_ice * g * thickness + rho_water * g * base, r = q / p and s = 1 / p)\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'coefficient', 'friction coefficient [SI]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'pressure_adjusted_temperature', 'friction pressure_adjusted_temperature (T - Tpmp) [K]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'gamma', '(T - Tpmp)/gamma [K]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.josh Class'
        return s
    
    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude friction.josh fields to 3D
        """
        self.coefficient = mesh.project_3d(md, vector = self.coefficient, type = 'node', layer = 1)
        self.pressure_adjusted_temperature = mesh.project_3d(md, vector = self.pressure_adjusted_temperature, type = 'node', layer = 1)
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [friction.josh] parameters.

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

        # Early return if necessary analyses or solutions not specified
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        
        class_utils.check_field(md, fieldname = "friction.coefficient", timeseries = True, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "friction.pressure_adjusted_temperature", allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "friction.gamma", gt = 0, scalar = True, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "friction.effective_pressure_limit", scalar = True, ge = 0)       

        class_utils.check_field(md, fieldname = "initialization.temperature", size = 'universal', allow_nan = False, allow_inf = False)

        return md
    
    # Marshall method for saving the friction.josh parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.josh] parameters to a binary file.

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

        ## Write headers to file
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.friction.law', data = 9, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'coefficient', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'pressure_adjusted_temperature', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'gamma', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'effective_pressure_limit', format = 'Double')

## ------------------------------------------------------
## friction.pism
## ------------------------------------------------------
@class_registry.register_class
class pism(class_registry.manage_state):
    """
    PISM friction class for ISSM.

    This class contains the parameters for the PISM friction law in the ISSM framework.
    It defines the main friction-related parameters specific to the PISM law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    pseudoplasticity_exponent : :class:`float`, default=0.6
        Pseudoplasticity exponent [dimensionless].
    threshold_speed : :class:`float`, default=100.
        Threshold speed [m/yr].
    delta : :class:`float`, default=0.02
        Lower limit of the effective pressure, expressed as a fraction of overburden pressure [dimensionless].
    void_ratio : :class:`float`, default=0.69
        Void ratio at a reference effective pressure [dimensionless].
    till_friction_angle : :class:`float`, default=np.nan
        Till friction angle [deg], recommended default: 30 deg.
    sediment_compressibility_coefficient : :class:`float`, default=np.nan
        Coefficient of compressibility of the sediment [dimensionless], recommended default: 0.12.

    Examples
    --------
    .. code-block:: python

        >>> md.friction = pyissm.model.classes.friction.pism()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.pseudoplasticity_exponent = 0.6
        self.threshold_speed = 100.
        self.delta = 0.02
        self.void_ratio = 0.69
        self.till_friction_angle = np.nan
        self.sediment_compressibility_coefficient = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Basal shear stress parameters for the PISM friction law (See Aschwanden et al. 2016 for more details)\n'

        s += "{}\n".format(class_utils.fielddisplay(self, 'pseudoplasticity_exponent', 'pseudoplasticity exponent [dimensionless]'))
        s += "{}\n".format(class_utils.fielddisplay(self, 'threshold_speed', 'threshold speed [m / yr]'))
        s += "{}\n".format(class_utils.fielddisplay(self, 'delta', 'lower limit of the effective pressure, expressed as a fraction of overburden pressure [dimensionless]'))
        s += "{}\n".format(class_utils.fielddisplay(self, 'void_ratio', 'void ratio at a reference effective pressure [dimensionless]'))
        s += "{}\n".format(class_utils.fielddisplay(self, 'till_friction_angle', 'till friction angle [deg], recommended default: 30 deg'))
        s += "{}\n".format(class_utils.fielddisplay(self, 'sediment_compressibility_coefficient', 'coefficient of compressibility of the sediment [dimensionless], recommended default: 0.12'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.pism Class'
        return s
    
    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude friction.pism fields to 3D
        """
        self.till_friction_angle = mesh.project_3d(md, vector = self.till_friction_angle, type = 'node', layer = 1)
        self.sediment_compressibility_coefficient = mesh.project_3d(md, vector = self.sediment_compressibility_coefficient, type = 'node', layer = 1)
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [friction.pism] parameters.

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

        # Early return if necessary analyses or solutions not specified
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        if solution == 'TransientSolution' and not md.transient.isstressbalance and not md.transient.isthermal:
            return md

        class_utils.check_field(md, fieldname = 'friction.pseudoplasticity_exponent', scalar = True, gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'friction.threshold_speed', scalar = True, gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'friction.delta', scalar = True, gt = 0, lt = 1, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'friction.void_ratio', scalar = True, gt = 0, lt = 1, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'friction.till_friction_angle', gt = 0, lt = 360., size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'friction.sediment_compressibility_coefficient', gt = 0., lt = 1., size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)

        return md

    # Marshall method for saving the friction.pism parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.pism] parameters to a binary file.

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

        ## Write headers to file
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.friction.law', data = 10, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'pseudoplasticity_exponent', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'threshold_speed', format = 'Double', scale = 1. / md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'delta', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'void_ratio', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'till_friction_angle', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'sediment_compressibility_coefficient', format = 'DoubleMat', mattype = 1)

## ------------------------------------------------------
## friction.regcoulomb
## ------------------------------------------------------
@class_registry.register_class
class regcoulomb(class_registry.manage_state):
    """
    Regularized Coulomb friction class for ISSM.

    This class contains the parameters for the regularized Coulomb friction law (Joughin et al., 2019) in the ISSM framework.
    It defines the main friction-related parameters specific to the regularized Coulomb law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    C : :class:`float` or :class:`numpy.ndarray`, default=np.nan
        Friction coefficient [SI].
    u0 : :class:`float` or :class:`numpy.ndarray`, default=1000
        Velocity controlling plastic limit.
    m : :class:`float` or :class:`numpy.ndarray`, default=np.nan
        m exponent (set to m = 3 in original paper).

    Examples
    --------
    .. code-block:: python
    
        >>> md.friction = pyissm.model.classes.friction.regcoulomb()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.C = np.nan
        self.u0 = 1000
        self.m = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        # See Joughin et al. 2019 (equivalent form by Matt Trevers, poster at AGU 2022) https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019GL082526
        s = 'Regularized Coulomb friction law (Joughin et al., 2019) parameters:\n'

        s += '   Regularized Coulomb friction law reads:\n'
        s += '                       C^2 |u|^(1/m)         \n'
        s += '      tau_b = -  ____________________________\n'
        s += '                     (|u|/u0 + 1)^(1/m)      \n'
        s += '\n'
        s += '{}\n'.format(class_utils.fielddisplay(self, 'C', 'friction coefficient [SI]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'm', 'm exponent (set to m = 3 in original paper)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'u0', 'velocity controlling plastic limit'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.regcoulomb Class'
        return s
    
    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude friction.regcoloumb fields to 3D
        """
        self.C = mesh.project_3d(md, vector = self.C, type = 'node')
        self.m = mesh.project_3d(md, vector = self.m, type = 'element')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [friction.regcoulomb] parameters.

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

        # Early return if necessary analyses or solutions not specified
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        
        class_utils.check_field(md, fieldname = 'friction.C', timeseries = True, ge = 0., allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'friction.u0', scalar = True, gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'friction.m', size = (md.mesh.numberofelements, 1), gt = 0., allow_nan = False, allow_inf = False)

        return md

    # Marshall method for saving the friction.regcoulomb parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.regcoulomb] parameters to a binary file.

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

        ## Write headers to file
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.friction.law', data = 14, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'C', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'u0', format = 'Double', scale = 1. / md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'm', format = 'DoubleMat', mattype = 2)

## ------------------------------------------------------
## friction.regcoulomb2
## ------------------------------------------------------
@class_registry.register_class
class regcoulomb2(class_registry.manage_state):
    """
    Regularized Coulomb 2 friction class for ISSM.

    This class contains the parameters for the regularized Coulomb 2 friction law (see Zoet and Iverson 2020 or Choi et al., 2022) in the ISSM framework.
    It defines the main friction-related parameters specific to the regularized Coulomb 2 law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    C : :class:`numpy.ndarray`, default=np.nan
        Friction coefficient [SI].
    K : :class:`numpy.ndarray`, default=np.nan
        K parameter for velocity controlling plastic limit.
    m : :class:`numpy.ndarray`, default=np.nan
        m exponent.
    effective_pressure_limit : :class:`float`, default=0
        Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0).

    Examples
    --------
    .. code-block:: python
        
        >>> md.friction = pyissm.model.classes.friction.regcoulomb2()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.C = np.nan
        self.K = np.nan
        self.m = np.nan
        self.effective_pressure_limit = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        # See Zoet and Iverson 2020 or Choi et al., 2022
        s = 'Regularized Coulomb friction law 2 parameters:\n'

        s += '   Regularized Coulomb friction law reads:\n'
        s += '                       C N |u|^(1/m)         \n'
        s += '      tau_b = -  ____________________________\n'
        s += '                   (|u| + (K*N)^m)^(1/m)     \n'
        s += '\n'
        s += '{}\n'.format(class_utils.fielddisplay(self, 'C', 'friction coefficient [SI]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'm', 'm exponent'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'K', '(K * N) ^ m to be velocity controlling plastic limit'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.regcoulomb2 Class'
        return s
    
    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude friction.regcoloumb2 fields to 3D
        """
        self.C = mesh.project_3d(md, vector = self.C, type = 'node')
        self.m = mesh.project_3d(md, vector = self.m, type = 'element')
        self.K = mesh.project_3d(md, vector = self.K, type = 'node')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [friction.regcoulomb2] parameters.

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

        # Early return if necessary analyses or solutions not specified
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        
        class_utils.check_field(md, fieldname = 'friction.C', timeseries = True, ge = 0., allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'friction.K', gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'friction.m', size = (md.mesh.numberofelements, 1), gt = 0., allow_nan = False, allow_inf = False)

        return md

    # Marshall method for saving the friction.regcoulomb2 parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.regcoulomb2] parameters to a binary file.

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

        ## Write headers to file
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.friction.law', data = 15, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'C', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'K', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'm', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'effective_pressure_limit', format = 'Double')

## ------------------------------------------------------
## friction.schoof
## ------------------------------------------------------
@class_registry.register_class
class schoof(class_registry.manage_state):
    """
    Schoof friction class for ISSM.

    This class contains the parameters for the Schoof sliding law in the ISSM framework.
    It defines the main friction-related parameters specific to the Schoof law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    C : :class:`numpy.ndarray`, default=np.nan
        Friction coefficient [SI].
    Cmax : :class:`numpy.ndarray`, default=np.nan
        Iken's bound (typically between 0.17 and 0.84) [SI].
    m : :class:`numpy.ndarray`, default=np.nan
        m exponent (generally taken as m = 1/n = 1/3).
    coupling : :class:`int`, default=0
        Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: used coupled model (not implemented yet).
    effective_pressure : :class:`numpy.ndarray`, default=np.nan
        Effective Pressure for the forcing if not coupled [Pa].
    effective_pressure_limit : :class:`float`, default=0
        Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0).

    Examples
    --------
    .. code-block:: python

        >>> md.friction = pyissm.model.classes.friction.schoof()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.C = np.nan
        self.Cmax = np.nan
        self.m = np.nan
        self.coupling = 0
        self.effective_pressure = np.nan
        self.effective_pressure_limit = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        # See Brondex et al. 2019 https://www.the-cryosphere.net/13/177/2019/
        s = 'Schoof sliding law parameters:\n'

        s += '   Schoof\'s sliding law reads:\n'
        s += '                         C^2 |u_b|^(m-1)                \n'
        s += '      tau_b = - _____________________________   u_b   \n'
        s += '               (1+(C^2/(Cmax N))^1/m |u_b| )^m          \n'
        s += '\n'
        s += "{}\n".format(class_utils.fielddisplay(self, 'C', 'friction coefficient [SI]'))
        s += "{}\n".format(class_utils.fielddisplay(self, 'Cmax', 'Iken\'s bound (typically between 0.17 and 0.84) [SI]'))
        s += "{}\n".format(class_utils.fielddisplay(self, 'm', 'm exponent (generally taken as m = 1/n = 1/3)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'coupling', 'Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: used coupled model (not implemented yet)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        s += "{}\n".format(class_utils.fielddisplay(self, 'effective_pressure_limit', 'fNeff do not allow to fall below a certain limit: effective_pressure_limit*rho_ice*g*thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.schoof Class'
        return s
    
    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude friction.schoof fields to 3D
        """
        self.C = mesh.project_3d(md, vector = self.C, type = 'node')
        self.Cmax = mesh.project_3d(md, vector = self.Cmax, type = 'node')
        self.m = mesh.project_3d(md, vector = self.m, type = 'element')
        if self.coupling in [3, 4]:
            self.effective_pressure = mesh.project_3d(md, vector = self.effective_pressure, type = 'node', layer = 1)
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [friction.schoof] parameters.

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

        # Early return if necessary analyses or solutions not specified
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        
        class_utils.check_field(md, fieldname = 'friction.C', timeseries = True, gt = 0., allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'friction.Cmax', timeseries = True, gt = 0., allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'friction.m', size = (md.mesh.numberofelements, 1), gt = 0., allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'friction.effective_pressure_limit', scalar = True, ge = 0.)
        class_utils.check_field(md, fieldname = 'friction.coupling', scalar = True, values = [0, 1, 2, 3, 4])
        
        if self.coupling == 3:
            class_utils.check_field(md, fieldname = 'friction.effective_pressure', timeseries = True, allow_nan = False, allow_inf = False)

        return md
    
    # Marshall method for saving the friction.schoof parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.schoof] parameters to a binary file.

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

        ## Write headers to file
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.friction.law', data = 11, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'C', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'Cmax', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'm', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'coupling', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'effective_pressure_limit', format = 'Double')

        ## Write conditional effective pressure
        if self.coupling in [3, 4]:
            execute.WriteData(fid, prefix, obj = self, fieldname = 'effective_pressure', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        elif self.coupling > 4:
            raise ValueError(f'md.friction.coupling = {self.coupling} is not implemented yet')

## ------------------------------------------------------
## friction.shakti
## ------------------------------------------------------
@class_registry.register_class
class shakti(class_registry.manage_state):
    """
    Shakti friction class for ISSM.

    This class contains the parameters for the Shakti friction law in the ISSM framework.
    It defines the main friction-related parameters specific to the Shakti law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    coefficient : :class:`numpy.ndarray`, default=np.nan
        Friction coefficient [SI].

    Examples
    --------
    .. code-block:: python
    
        >>> md.friction = pyissm.model.classes.friction.shakti()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.coefficient = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Basal shear stress parameters: Sigma_b = coefficient^2 * Neff * u_b\n'

        s += '(effective stress Neff = rho_ice * g * thickness + rho_water * g * (head - b))\n'
        s += '{}\n'.format(class_utils.fielddisplay(self, 'coefficient', 'friction coefficient [SI]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.shakti Class'
        return s
    
    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude friction.shakti fields to 3D
        """
        self.coefficient = mesh.project_3d(md, vector = self.coefficient, type = 'node', layer = 1)
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [friction.shakti] parameters.

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

        # Early return if necessary analyses or solutions not specified
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        
        class_utils.check_field(md, fieldname = "friction.coefficient", timeseries = True, allow_nan = False, allow_inf = False)

        return md
    
    # Marshall method for saving the friction.shakti parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.shakti] parameters to a binary file.

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

        ## Write headers to file
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.friction.law', data = 8, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'coefficient', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

## ------------------------------------------------------
## friction.waterlayer
## ------------------------------------------------------
@class_registry.register_class
class waterlayer(class_registry.manage_state):
    """
    Waterlayer friction class for ISSM.

    This class contains the parameters for the waterlayer friction law in the ISSM framework.
    It defines the main friction-related parameters specific to the waterlayer law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    coefficient : :class:`numpy.ndarray`, default=np.nan
        Friction coefficient [SI].
    f : :class:`numpy.ndarray`, default=np.nan
        f variable for effective pressure.
    p : :class:`numpy.ndarray`, default=np.nan
        p exponent.
    q : :class:`numpy.ndarray`, default=np.nan
        q exponent.
    water_layer : :class:`numpy.ndarray`, default=np.nan
        Water thickness at the base of the ice (m).

    Examples
    --------
    .. code-block:: python
    
        >>> md.friction = pyissm.model.classes.friction.waterlayer()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.coefficient = np.nan
        self.f = np.nan
        self.p = np.nan
        self.q = np.nan
        self.water_layer = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Basal shear stress parameters: tau_b = coefficient^2 * Neff ^r * |u_b|^(s - 1) * u_b * 1 / f(T)\n(effective stress Neff = rho_ice * g * thickness + rho_water * g * (bed + water_layer), r = q / p and s = 1 / p)\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'coefficient', 'frictiontemp coefficient [SI]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'f', 'f variable for effective pressure'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'p', 'p exponent'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'q', 'q exponent'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'water_layer', 'water thickness at the base of the ice (m)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.waterlayer Class'
        return s

    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude friction.schoof fields to 3D
        """
        self.coefficient = mesh.project_3d(md, vector = self.coefficient, type = 'node', layer = '1')
        self.p = mesh.project_3d(md, vector = self.p, type = 'element')
        self.p = mesh.project_3d(md, vector = self.q, type = 'element')
        self.water_layer = mesh.project_3d(md, vector = self.water_layer, type = 'node', layer = 1)
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [friction.waterlayer] parameters.

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

        # Early return if necessary analyses or solutions not specified
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        
        class_utils.check_field(md, fieldname = 'friction.coefficient', timeseries = True, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'friction.f', scalar = True, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'friction.q', size = (md.mesh.numberofelements, ), allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'friction.p', size = (md.mesh.numberofelements, ), allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'thermal.spctemperature', timeseries = True, ge = 0., allow_inf = False)

        return md

    # Marshall method for saving the friction.waterlayer parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.waterlayer] parameters to a binary file.

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

        ## Write headers to file
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.friction.law', data = 5, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'coefficient', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'f', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'p', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'q', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'water_layer', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        
## ------------------------------------------------------
## friction.weertman
## ------------------------------------------------------
@class_registry.register_class
class weertman(class_registry.manage_state):
    """
    Weertman friction class for ISSM.

    This class contains the parameters for the Weertman sliding law in the ISSM framework.
    It defines the main friction-related parameters specific to the Weertman law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    C : :class:`numpy.ndarray`, default=np.nan
        Friction coefficient [SI].
    m : :class:`numpy.ndarray`, default=np.nan
        m exponent.
    linearize : :class:`int`, default=0
        0: not linearized, 1: interpolated linearly, 2: constant per element (default is 0).

    Examples
    --------
    .. code-block:: python
    
        >>> md.friction = pyissm.model.classes.friction.weertman()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.C = np.nan
        self.m = np.nan
        self.linearize = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Weertman sliding law parameters: Sigma_b = C^(- 1 / m) * |u_b|^(1 / m - 1) * u_b\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'C', 'friction coefficient [SI]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'm', 'm exponent'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'linearize', '0: not linearized, 1: interpolated linearly, 2: constant per element (default is 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.weertman Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [friction.weertman] parameters.

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

        # Early return if necessary analyses or solutions not specified
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        
        class_utils.check_field(md, fieldname = 'friction.C', timeseries = True, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'friction.m', size = (md.mesh.numberofelements, ), allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'friction.linearize', scalar = True, values = [0, 1, 2])

        return md

    # Marshall method for saving the friction.weertman parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.weertman] parameters to a binary file.

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

        ## Write headers to file
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.friction.law', data = 2, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'C', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'm', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'linearize', format = 'Integer')
