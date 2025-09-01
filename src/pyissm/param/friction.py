import numpy as np
from . import param_utils
from . import class_registry
from .. import execute

## ------------------------------------------------------
## friction.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    """
    Default friction parameters class for ISSM.

    This class encapsulates the default parameters for friction in the ISSM (Ice Sheet System Model) framework.
    It defines the main friction-related parameters.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    coefficient : ndarray, default=np.nan
        Friction coefficient [SI].
    p : ndarray, default=np.nan
        p exponent.
    q : ndarray, default=np.nan
        q exponent.
    coupling : int, default=0
        Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: used coupled model (not implemented yet).
    linearize : int, default=0
        0: not linearized, 1: interpolated linearly, 2: constant per element (default is 0).
    effective_pressure : ndarray, default=np.nan
        Effective Pressure for the forcing if not coupled [Pa].
    effective_pressure_limit : ndarray, default=0
        Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0).

    Methods
    -------
    __init__(self, other=None)
        Initializes the friction parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the friction parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.friction = pyissm.param.friction.default()
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

        s += '{}\n'.format(param_utils.fielddisplay(self, 'coefficient', 'friction coefficient [SI]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'p', 'p exponent'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'q', 'q exponent'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'coupling', 'Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: used coupled model (not implemented yet)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'linearize', '0: not linearized, 1: interpolated linearly, 2: constant per element (default is 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.default Class'
        return s


    # Marshall method for saving the friction.default parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.default] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
            ISSM model object needed in some cases.

        Returns
        -------
        None
        """

        ## Write headers to file
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, data = 1, name = 'md.friction.law', format = 'Integer')

        ## Write coefficient field
        if isinstance(self.coefficient, np.ndarray) and ((self.coefficient.shape[0] == md.mesh.numberofvertices) or self.coefficient.shape[0] == md.mesh.numberofvertices + 1):
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
    Coulomb friction parameters class for ISSM.

    This class encapsulates the parameters for the Coulomb friction law in the ISSM (Ice Sheet System Model) framework.
    It defines the main friction-related parameters specific to the Coulomb law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    coefficient : ndarray, default=np.nan
        Power law (Weertman) friction coefficient [SI].
    coefficientcoulomb : ndarray, default=np.nan
        Coulomb friction coefficient [SI].
    p : ndarray, default=np.nan
        p exponent.
    q : ndarray, default=np.nan
        q exponent.
    coupling : int, default=0
        Coupling flag: 0 for default, 1 for forcing (provide md.friction.effective_pressure), 2 for coupled (not implemented yet).
    effective_pressure : =ndarray, default=np.nan
        Effective Pressure for the forcing if not coupled [Pa].
    effective_pressure_limit : =ndarray, default=0
        Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0).

    Methods
    -------
    __init__(self, other=None)
        Initializes the friction parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the friction parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.friction = pyissm.param.friction.coulomb()
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

        s += '{}\n'.format(param_utils.fielddisplay(self, 'coefficient', 'power law (Weertman) friction coefficient [SI]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'coefficientcoulomb', 'Coulomb friction coefficient [SI]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'p', 'p exponent'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'q', 'q exponent'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'coupling', 'Coupling flag: 0 for default, 1 for forcing(provide md.friction.effective_pressure)  and 2 for coupled(not implemented yet)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.coulomb Class'
        return s
    
    # Marshall method for saving the friction.coulomb parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.coulomb] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
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
    Coulomb2 friction parameters class for ISSM.

    This class encapsulates the parameters for the Coulomb2 friction law in the ISSM (Ice Sheet System Model) framework.
    It defines the main friction-related parameters specific to the Coulomb2 law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    coefficient : ndarray, default=np.nan
        Power law (Weertman) friction coefficient [SI].
    coefficientcoulomb : ndarray, default=np.nan
        Coulomb friction coefficient [SI].
    p : ndarray, default=np.nan
        p exponent.
    q : ndarray, default=np.nan
        q exponent.
    coupling : int, default=0
        Coupling flag: 0 for default, 1 for forcing (provide md.friction.effective_pressure), 2 for coupled (not implemented yet).
    effective_pressure : ndarray, default=np.nan
        Effective Pressure for the forcing if not coupled [Pa].
    effective_pressure_limit : ndarray, default=0
        Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0).

    Methods
    -------
    __init__(self, other=None)
        Initializes the friction parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the friction parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.friction = pyissm.param.friction.coulomb2()
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

        s += '{}\n'.format(param_utils.fielddisplay(self, 'coefficient', 'power law (Weertman) friction coefficient [SI]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'coefficientcoulomb', 'Coulomb friction coefficient [SI]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'p', 'p exponent'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'q', 'q exponent'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'coupling', 'Coupling flag: 0 for default, 1 for forcing(provide md.friction.effective_pressure)  and 2 for coupled(not implemented yet)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.coulomb2 Class'
        return s

    # Marshall method for saving the friction.coulomb2 parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.coulomb2] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
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
    Hydro friction parameters class for ISSM.

    This class encapsulates the parameters for the hydro (Gagliardini 2007) friction law in the ISSM (Ice Sheet System Model) framework.
    It defines the main friction-related parameters specific to the hydro law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    coupling : int, default=0
        Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: use coupled model (not implemented yet).
    q : ndarray, default=np.nan
        Friction law exponent q >= 1.
    C : ndarray, default=np.nan
        Friction law max value (Iken bound).
    As : ndarray, default=np.nan
        Sliding parameter without cavitation [m Pa^-n s^-1].
    effective_pressure : ndarray, default=np.nan
        Effective Pressure for the forcing if not coupled [Pa].
    effective_pressure_limit : ndarray, default=0
        Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0).

    Methods
    -------
    __init__(self, other=None)
        Initializes the friction parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the friction parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.friction = pyissm.param.friction.hydro()
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

        s += '{}\n'.format(param_utils.fielddisplay(self, 'coupling', 'Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: used coupled model (not implemented yet)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'q', 'friction law exponent q >= 1'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'C', 'friction law max value (Iken bound)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'As', 'Sliding Parameter without cavitation [m Pa^ - n s^ - 1]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.hydro Class'
        return s

    # Marshall method for saving the friction.hydro parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.hydro] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
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
    Josh friction parameters class for ISSM.

    This class encapsulates the parameters for the Josh friction law in the ISSM (Ice Sheet System Model) framework.
    It defines the main friction-related parameters specific to the Josh law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    coefficient : ndarray, default=np.nan
        Friction coefficient [SI].
    pressure_adjusted_temperature : ndarray, default=np.nan
        Friction pressure_adjusted_temperature (T - Tpmp) [K].
    gamma : ndarray, default=1.
        (T - Tpmp)/gamma [K].
    effective_pressure_limit : ndarray, default=0
        Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0).

    Methods
    -------
    __init__(self, other=None)
        Initializes the friction parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the friction parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.friction = pyissm.param.friction.josh()
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

        s += '{}\n'.format(param_utils.fielddisplay(self, 'coefficient', 'friction coefficient [SI]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'pressure_adjusted_temperature', 'friction pressure_adjusted_temperature (T - Tpmp) [K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'gamma', '(T - Tpmp)/gamma [K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.josh Class'
        return s
    
    # Marshall method for saving the friction.josh parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.josh] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
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
    PISM friction parameters class for ISSM.

    This class encapsulates the parameters for the PISM friction law in the ISSM (Ice Sheet System Model) framework.
    It defines the main friction-related parameters specific to the PISM law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    pseudoplasticity_exponent : float, default=0.6
        Pseudoplasticity exponent [dimensionless].
    threshold_speed : float, default=100.
        Threshold speed [m/yr].
    delta : float, default=0.02
        Lower limit of the effective pressure, expressed as a fraction of overburden pressure [dimensionless].
    void_ratio : float, default=0.69
        Void ratio at a reference effective pressure [dimensionless].
    till_friction_angle : float, default=np.nan
        Till friction angle [deg], recommended default: 30 deg.
    sediment_compressibility_coefficient : float, default=np.nan
        Coefficient of compressibility of the sediment [dimensionless], recommended default: 0.12.
    requested_outputs : str, default='List of requested outputs'
        List of requested outputs.

    Methods
    -------
    __init__(self, other=None)
        Initializes the friction parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the friction parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.friction = pyissm.param.friction.pism()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.pseudoplasticity_exponent = 0.6
        self.threshold_speed = 100.
        self.delta = 0.02
        self.void_ratio = 0.69
        self.till_friction_angle = np.nan
        self.sediment_compressibility_coefficient = np.nan
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = 'Basal shear stress parameters for the PISM friction law (See Aschwanden et al. 2016 for more details)\n'

        s += "{}\n".format(param_utils.fielddisplay(self, 'pseudoplasticity_exponent', 'pseudoplasticity exponent [dimensionless]'))
        s += "{}\n".format(param_utils.fielddisplay(self, 'threshold_speed', 'threshold speed [m / yr]'))
        s += "{}\n".format(param_utils.fielddisplay(self, 'delta', 'lower limit of the effective pressure, expressed as a fraction of overburden pressure [dimensionless]'))
        s += "{}\n".format(param_utils.fielddisplay(self, 'void_ratio', 'void ratio at a reference effective pressure [dimensionless]'))
        s += "{}\n".format(param_utils.fielddisplay(self, 'till_friction_angle', 'till friction angle [deg], recommended default: 30 deg'))
        s += "{}\n".format(param_utils.fielddisplay(self, 'sediment_compressibility_coefficient', 'coefficient of compressibility of the sediment [dimensionless], recommended default: 0.12'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.pism Class'
        return s

    # Marshall method for saving the friction.pism parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.pism] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
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
    Regularized Coulomb friction parameters class for ISSM.

    This class encapsulates the parameters for the regularized Coulomb friction law (Joughin et al., 2019) in the ISSM (Ice Sheet System Model) framework.
    It defines the main friction-related parameters specific to the regularized Coulomb law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    C : float or ndarray, default=np.nan
        Friction coefficient [SI].
    u0 : float or ndarray, default=1000
        Velocity controlling plastic limit.
    m : float or ndarray, default=np.nan
        m exponent (set to m = 3 in original paper).

    Methods
    -------
    __init__(self, other=None)
        Initializes the friction parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the friction parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.friction = pyissm.param.friction.regcoulomb()
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
        s += '{}\n'.format(param_utils.fielddisplay(self, 'C', 'friction coefficient [SI]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'm', 'm exponent (set to m = 3 in original paper)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'u0', 'velocity controlling plastic limit'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.regcoulomb Class'
        return s
    
    # Marshall method for saving the friction.regcoulomb parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.regcoulomb] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
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
    Regularized Coulomb 2 friction parameters class for ISSM.

    This class encapsulates the parameters for the regularized Coulomb 2 friction law (see Zoet and Iverson 2020 or Choi et al., 2022) in the ISSM (Ice Sheet System Model) framework.
    It defines the main friction-related parameters specific to the regularized Coulomb 2 law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    C : \ndarray, default=np.nan
        Friction coefficient [SI].
    K : ndarray, default=np.nan
        K parameter for velocity controlling plastic limit.
    m : ndarray, default=np.nan
        m exponent.
    effective_pressure_limit : float, default=0
        Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0).

    Methods
    -------
    __init__(self, other=None)
        Initializes the friction parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the friction parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.friction = pyissm.param.friction.regcoulomb2()
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
        s += '{}\n'.format(param_utils.fielddisplay(self, 'C', 'friction coefficient [SI]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'm', 'm exponent'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'K', '(K * N) ^ m to be velocity controlling plastic limit'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.regcoulomb2 Class'
        return s

    # Marshall method for saving the friction.regcoulomb2 parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.regcoulomb2] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
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
    Schoof friction parameters class for ISSM.

    This class encapsulates the parameters for the Schoof sliding law in the ISSM (Ice Sheet System Model) framework.
    It defines the main friction-related parameters specific to the Schoof law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    C : ndarray, default=np.nan
        Friction coefficient [SI].
    Cmax : ndarray, default=np.nan
        Iken's bound (typically between 0.17 and 0.84) [SI].
    m : ndarray, default=np.nan
        m exponent (generally taken as m = 1/n = 1/3).
    coupling : int, default=0
        Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: used coupled model (not implemented yet).
    effective_pressure : ndarray, default=np.nan
        Effective Pressure for the forcing if not coupled [Pa].
    effective_pressure_limit : float, default=0
        Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0).

    Methods
    -------
    __init__(self, other=None)
        Initializes the friction parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the friction parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.friction = pyissm.param.friction.schoof()
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
        s += "{}\n".format(param_utils.fielddisplay(self, 'C', 'friction coefficient [SI]'))
        s += "{}\n".format(param_utils.fielddisplay(self, 'Cmax', 'Iken\'s bound (typically between 0.17 and 0.84) [SI]'))
        s += "{}\n".format(param_utils.fielddisplay(self, 'm', 'm exponent (generally taken as m = 1/n = 1/3)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'coupling', 'Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: used coupled model (not implemented yet)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        s += "{}\n".format(param_utils.fielddisplay(self, 'effective_pressure_limit', 'fNeff do not allow to fall below a certain limit: effective_pressure_limit*rho_ice*g*thickness (default 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.schoof Class'
        return s
    
    # Marshall method for saving the friction.schoof parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.schoof] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
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
    Shakti friction parameters class for ISSM.

    This class encapsulates the parameters for the Shakti friction law in the ISSM (Ice Sheet System Model) framework.
    It defines the main friction-related parameters specific to the Shakti law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    coefficient : ndarray, default=np.nan
        Friction coefficient [SI].

    Methods
    -------
    __init__(self, other=None)
        Initializes the friction parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the friction parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.friction = pyissm.param.friction.shakti()
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
        s += '{}\n'.format(param_utils.fielddisplay(self, 'coefficient', 'friction coefficient [SI]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.shakti Class'
        return s
    
    # Marshall method for saving the friction.shakti parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.shakti] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
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
    Waterlayer friction parameters class for ISSM.

    This class encapsulates the parameters for the waterlayer friction law in the ISSM (Ice Sheet System Model) framework.
    It defines the main friction-related parameters specific to the waterlayer law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    coefficient : ndarray, default=np.nan
        Friction coefficient [SI].
    f : ndarray, default=np.nan
        f variable for effective pressure.
    p : ndarray, default=np.nan
        p exponent.
    q : ndarray, default=np.nan
        q exponent.
    water_layer : ndarray, default=np.nan
        Water thickness at the base of the ice (m).

    Methods
    -------
    __init__(self, other=None)
        Initializes the friction parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the friction parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.friction = pyissm.param.friction.waterlayer()
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

        s += '{}\n'.format(param_utils.fielddisplay(self, 'coefficient', 'frictiontemp coefficient [SI]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'f', 'f variable for effective pressure'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'p', 'p exponent'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'q', 'q exponent'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'water_layer', 'water thickness at the base of the ice (m)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.waterlayer Class'
        return s

    # Marshall method for saving the friction.waterlayer parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.waterlayer] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
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
    Weertman friction parameters class for ISSM.

    This class encapsulates the parameters for the Weertman sliding law in the ISSM (Ice Sheet System Model) framework.
    It defines the main friction-related parameters specific to the Weertman law.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    C : ndarray, default=np.nan
        Friction coefficient [SI].
    m : ndarray, default=np.nan
        m exponent.
    linearize : int, default=0
        0: not linearized, 1: interpolated linearly, 2: constant per element (default is 0).

    Methods
    -------
    __init__(self, other=None)
        Initializes the friction parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the friction parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.friction = pyissm.param.friction.weertman()
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

        s += '{}\n'.format(param_utils.fielddisplay(self, 'C', 'friction coefficient [SI]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'm', 'm exponent'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'linearize', '0: not linearized, 1: interpolated linearly, 2: constant per element (default is 0)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - friction.weertman Class'
        return s

    # Marshall method for saving the friction.weertman parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [friction.weertman] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
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
