import numpy as np
import warnings
from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm.model import execute

## ------------------------------------------------------
## basalforcings.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    """
    Default basal forcings parameters class for ISSM.

    This class encapsulates the default parameters for basal forcings in the ISSM (Ice Sheet System Model) framework.
    It defines the melting rates for grounded and floating ice, perturbation melting rate, and geothermal heat flux.
    
    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    groundedice_melting_rate : ndarray, default=np.nan
        Basal melting rate for grounded ice (positive if melting) [m/yr].
    floatingice_melting_rate : ndarray, default=np.nan
        Basal melting rate for floating ice (positive if melting) [m/yr].
    perturbation_melting_rate : ndarray, default=np.nan
        Optional perturbation in basal melting rate under floating ice (positive if melting) [m/yr].
    geothermalflux : float, default=np.nan
        Geothermal heat flux [W/m^2].

    Methods
    -------
    __init__(self, other=None)
        Initializes the basal forcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the basal forcings parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.basalforcings = pyissm.model.classes.basalforcings.default()
    md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices,))
    md.basalforcings.floatingice_melting_rate = np.ones((md.mesh.numberofvertices,)) * 2
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.groundedice_melting_rate = np.nan
        self.floatingice_melting_rate = np.nan
        self.perturbation_melting_rate = np.nan
        self.geothermalflux = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   basal forcings parameters:\n'
        
        s += '{}\n'.format(class_utils.fielddisplay(self, 'groundedice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'floatingice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'perturbation_melting_rate', '(optional) perturbation in basal melting rate under floating ice [m/yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'geothermalflux', 'geothermal heat flux [W/m^2]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - basalforcings.default Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        # Masstransport analysis
        if 'Masstransport' in analyses and solution != 'TransientSolution' and not md.transient.ismasstransport:
            class_utils.check_field(md, fieldname = "basalforcings.groundedice_melting_rate", timeseries = True, allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = "basalforcings.floatingice_melting_rate", timeseries = True, allow_nan = False, allow_inf = False)

        # BalancethicknessAnalysis
        if 'BalancethicknessAnalysis' in analyses:
            class_utils.check_field(md, fieldname = "basalforcings.groundedice_melting_rate", size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = "basalforcings.floatingice_melting_rate", size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)

        if 'ThermalAnalysis' in analyses and solution != 'TransientSolution' and not md.transient.isthermal:
            class_utils.check_field(md, fieldname = "basalforcings.groundedice_melting_rate", timeseries = True, allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = "basalforcings.floatingice_melting_rate", timeseries = True, allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = "basalforcings.geothermalflux", timeseries = True, ge = 0, allow_nan = False, allow_inf = False)
            
        return md
    
    # Initialise empty fields of correct dimensions
    def initialise(self, md):
        """
        Initialise empty fields in basalforcings.default.
        """

        if np.all(np.isnan(self.groundedice_melting_rate)):
            self.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices,))
            warnings.warn('pyissm.model.classes.basalforcings.default: no groundedice_melting_rate specified -- values set as 0')

        if np.all(np.isnan(self.floatingice_melting_rate)):
            self.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices,))
            warnings.warn('pyissm.model.classes.basalforcings.default: no floatingice_melting_rate specified -- values set as 0')

        return self

    # Marshall method for saving the basalforcings.default parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [basalforcings.default] parameters to a binary file.

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

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.basalforcings.model', data = 1, format = 'Integer')

        ## Write scaled fields
        fieldnames = ['groundedice_melting_rate', 'floatingice_melting_rate', 'perturbation_melting_rate']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'geothermalflux', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

## ------------------------------------------------------
## basalforcings.pico
## ------------------------------------------------------
@class_registry.register_class
class pico(class_registry.manage_state):
    """
    Potsdam Ice-shelf Cavity mOdel (PICO) basal forcings parameters class for ISSM.

    This class encapsulates the parameters for the PICO basal melt parameterization in the ISSM (Ice Sheet System Model) framework.
    It defines the structure of the ice shelf cavities, including the number of basins, basin IDs, and various parameters related to ocean temperature, salinity, and melting rates.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    num_basins : int, default=0
        Number of basins the model domain is partitioned into [unitless].
    basin_id : ndarray, default=np.nan
        Basin number assigned to each node [unitless].
    maxboxcount : int, default=0
        Maximum number of boxes initialized under all ice shelves.
    overturning_coeff : float, default=np.nan
        Overturning strength [m^3/s].
    gamma_T : float, default=0.
        Turbulent temperature exchange velocity [m/s].
    farocean_temperature : ndarray, default=np.nan
        Depth averaged ocean temperature in front of the ice shelf for each basin [K].
    farocean_salinity : ndarray, default=np.nan
        Depth averaged ocean salinity in front of the ice shelf for each basin [psu].
    isplume : int, default=0
        Boolean (0 or 1) to use buoyant plume melt rate parameterization from Lazeroms et al., 2018 (default false).
    geothermalflux : float, default=np.nan
        Geothermal heat flux [W/m^2].
    groundedice_melting_rate : ndarray, default=np.nan
        Basal melting rate for grounded ice (positive if melting) [m/yr].

    Methods
    -------
    __init__(self, other=None)
        Initializes the PICO basal forcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the PICO basal forcings parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.basalforcings = pyissm.model.classes.basalforcings.pico()
    md.basalforcings.num_basins = 3
    md.basalforcings.basin_id = np.array([1, 2, 3])
    md.basalforcings.farocean_temperature = np.array([273.15, 273.2, 273.1])
    """

    # Initialise with default parameters
    def __init__(self, other=None):
        self.num_basins = 0
        self.basin_id = np.nan
        self.maxboxcount = 0
        self.overturning_coeff = np.nan
        self.gamma_T = 0.
        self.farocean_temperature = np.nan
        self.farocean_salinity = np.nan
        self.isplume = 0
        self.geothermalflux = np.nan
        self.groundedice_melting_rate = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   PICO basal melt rate parameterization:\n'

        s += '{}\n'.format(class_utils.fielddisplay(self,'num_basins','number of basins the model domain is partitioned into [unitless]'))
        s += '{}\n'.format(class_utils.fielddisplay(self,'basin_id','basin number assigned to each node [unitless]'))
        s += '{}\n'.format(class_utils.fielddisplay(self,'maxboxcount','maximum number of boxes initialized under all ice shelves'))
        s += '{}\n'.format(class_utils.fielddisplay(self,'overturning_coeff','overturning strength [m^3/s]'))
        s += '{}\n'.format(class_utils.fielddisplay(self,'gamma_T','turbulent temperature exchange velocity [m/s]'))
        s += '{}\n'.format(class_utils.fielddisplay(self,'farocean_temperature','depth averaged ocean temperature in front of the ice shelf for basin i [K]'))
        s += '{}\n'.format(class_utils.fielddisplay(self,'farocean_salinity','depth averaged ocean salinity in front of the ice shelf for basin i [psu]'))
        s += '{}\n'.format(class_utils.fielddisplay(self,'isplume','boolean to use buoyant plume melt rate parameterization from Lazeroms et al., 2018 (default false)'))
        s += '{}\n'.format(class_utils.fielddisplay(self,'geothermalflux','geothermal heat flux [W/m^2]'))
        s += '{}\n'.format(class_utils.fielddisplay(self,'groundedice_melting_rate','basal melting rate (positive if melting) [m/yr]'))

        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - basalforcings.pico Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):

        class_utils.check_field(md, fieldname = "basalforcings.num_basins", scalar = True, gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "basalforcings.basin_id", size = (md.mesh.numberofelements, 1), ge = 0, le = md.basalforcings.num_basins, allow_inf = False)
        class_utils.check_field(md, fieldname = "basalforcings.maxboxcount", scalar = True, gt = 0, allow_nan = False, allow_inf = False)

        if np.size(self.overturning_coeff) == 1:
            class_utils.check_field(md, fieldname = "basalforcings.overturning_coeff", scalar = True, gt = 0, allow_nan = False, allow_inf = False)
        else:
            class_utils.check_field(md, fieldname = "basalforcings.overturning_coeff", size = (md.mesh.numberofvertices, 1), gt = 0, allow_nan = False, allow_inf = False)

        class_utils.check_field(md, fieldname = "basalforcings.gamma_T", scalar = True, gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "basalforcings.farocean_temperature", size = (md.basalforcings.num_basins + 1, None), allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "basalforcings.farocean_salinity", size = (md.basalforcings.num_basins + 1, None), gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "basalforcings.isplume", scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = "basalforcings.geothermalflux", timeseries = True, ge = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "basalforcings.groundedice_melting_rate", timeseries = True, allow_nan = False, allow_inf = False)
    
        return md
    
    # Initialise empty fields of correct dimensions
    def initialise(self, md):
        """
        Initialise empty fields in basalforcings.pico.
        """

        if np.isnan(self.maxboxcount):
            self.maxboxcount = 5
            warnings.warn('pyissm.parm.basalforcings.pico: no maximum number of boxes set -- value set to 5.')

        if np.isnan(self.overturning_coeff):
            self.overturning_coeff = 1e6 * np.ones((md.mesh.numberofvertices,1))
            warnings.warn('pyissm.parm.basalforcings.pico: no overturning strength set -- value set to 1e6.')

        if np.isnan(self.gamma_T):
            self.gamma_T = 2e-5
            warnings.warn('pyissm.parm.basalforcings.pico: no turbulent temperature exchange velocity set -- value set to 2e-5.')

        if np.isnan(self.groundedice_melting_rate):
            self.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices,1))
            warnings.warn('pyissm.parm.basalforcings.pico: no basalforcings.groundedice_melting_rate specified -- values set as zero.')

        return self

    # Marshall method for saving the basalforcings.pico parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [basalforcings.pico] parameters to a binary file.

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

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.basalforcings.model', data = 5, format = 'Integer')

        ## Write Integer fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'num_basins', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'maxboxcount', format = 'Integer')

        ## Write Boolean fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'isplume', format = 'Boolean')

        ## Write DoubleMat fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'overturning_coeff', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'farocean_temperature', format = 'DoubleMat', timeserieslength = self.num_basins + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'farocean_salinity', format = 'DoubleMat', timeserieslength = self.num_basins + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'geothermalflux', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofelements + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'groundedice_melting_rate', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

        ## Write Double fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'gamma_T', format = 'Double')

        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'basin_id', data = self.basin_id - 1, format = 'IntMat', mattype = 2) # Change to 0-indexing

## ------------------------------------------------------
## basalforcings.linear
## ------------------------------------------------------
@class_registry.register_class
class linear(class_registry.manage_state):
    """
    Linear basal forcings parameters class for ISSM.

    This class encapsulates the parameters for linear basal forcings in the ISSM (Ice Sheet System Model) framework.
    It defines the melting rates for deep and upper water, grounded ice, and geothermal flux, allowing for a depth-dependent representation of basal melting processes.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    deepwater_melting_rate : float, default=0.
        Basal melting rate applied for floating ice with base < deepwater_elevation [m/yr].
    deepwater_elevation : float, default=0.
        Elevation threshold for deepwater melting rate [m].
    upperwater_melting_rate : float, default=0.
        Basal melting rate applied for floating ice with base >= upperwater_elevation [m/yr].
    upperwater_elevation : float, default=0.
        Elevation threshold for upperwater melting rate [m].
    groundedice_melting_rate : ndarray, default=np.nan
        Basal melting rate for grounded ice (positive if melting) [m/yr].
    perturbation_melting_rate : ndarray, default=np.nan
        Perturbation applied to computed melting rate (positive if melting) [m/yr].
    geothermalflux : float, default=np.nan
        Geothermal heat flux [W/m^2].

    Methods
    -------
    __init__(self, other=None)
        Initializes the linear basal forcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the linear basal forcings parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.basalforcings = pyissm.model.classes.basalforcings.linear()
    md.basalforcings.deepwater_melting_rate = 1.5
    md.basalforcings.deepwater_elevation = -500
    md.basalforcings.upperwater_melting_rate = 0.5
    md.basalforcings.upperwater_elevation = -200
    """

    # Initialise with default parameters
    def __init__(self, other=None):
        self.deepwater_melting_rate = 0.
        self.deepwater_elevation = 0.
        self.upperwater_melting_rate = 0.
        self.upperwater_elevation = 0.
        self.groundedice_melting_rate = np.nan
        self.perturbation_melting_rate = np.nan
        self.geothermalflux = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   linear basal forcings parameters:\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, "deepwater_melting_rate", "basal melting rate (positive if melting applied for floating ice whith base < deepwater_elevation) [m/yr]"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "deepwater_elevation", "elevation of ocean deepwater [m]"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "upperwater_melting_rate", "upper melting rate (positive if melting applied for floating ice whith base >= upperwater_elevation) [m/yr]"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "upperwater_elevation", "elevation of ocean upper water [m]"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "groundedice_melting_rate", "basal melting rate (positive if melting) [m/yr]"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "perturbation_melting_rate", "perturbation applied to computed melting rate (positive if melting) [m/yr]"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "geothermalflux", "geothermal heat flux [W/m^2]"))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - basalforcings.linear Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        if not np.all(np.isnan(self.perturbation_melting_rate)):
            class_utils.check_field(md, fieldname = "basalforcings.perturbation_melting_rate", timeseries = True, allow_nan = False, allow_inf = False)
        if 'MasstransportAnalysis' in analyses and solution != 'TransientSolution' and not md.transient.ismasstransport:
            class_utils.check_field(md, fieldname = "basalforcings.groundedice_melting_rate", timeseries = True, allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = "basalforcings.deepwater_melting_rate", singletimeseries = True, ge = 0)
            class_utils.check_field(md, fieldname = "basalforcings.upperwater_melting_rate", singletimeseries = True, ge = 0)
            class_utils.check_field(md, fieldname = "basalforcings.deepwater_elevation", singletimeseries = True, lt = self.upperwater_elevation)
            class_utils.check_field(md, fieldname = "basalforcings.upperwater_elevation", singletimeseries = True, le = 0)
        if 'BalancethicknessAnalysis' in analyses:
            class_utils.check_field(md, fieldname = "basalforcings.groundedice_melting_rate", size = (md.mesh.numberofvertices,), allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = "basalforcings.deepwater_melting_rate", singletimeseries = True, ge = 0)
            class_utils.check_field(md, fieldname = "basalforcings.upperwater_melting_rate", singletimeseries = True, ge = 0)
            class_utils.check_field(md, fieldname = "basalforcings.deepwater_elevation", singletimeseries = True, lt = self.upperwater_elevation)
            class_utils.check_field(md, fieldname = "basalforcings.upperwater_elevation", singletimeseries = True, le = 0)
        if 'ThermalAnalysis' in analyses and solution != 'TransientSolution' and not md.transient.isthermal:
            class_utils.check_field(md, fieldname = "basalforcings.groundedice_melting_rate", timeseries = True, allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = "basalforcings.deepwater_melting_rate", singletimeseries = True, ge = 0)
            class_utils.check_field(md, fieldname = "basalforcings.upperwater_melting_rate", singletimeseries = True, ge = 0)
            class_utils.check_field(md, fieldname = "basalforcings.deepwater_elevation", singletimeseries = True, lt = self.upperwater_elevation)
            class_utils.check_field(md, fieldname = "basalforcings.upperwater_elevation", singletimeseries = True, le = 0)
            class_utils.check_field(md, fieldname = "basalforcings.geothermalflux", timeseries = True, ge = 0, allow_nan = False, allow_inf = False)

        return md

    # Initialise empty fields of correct dimensions
    def initialise(self, md):
        """
        Initialise empty fields in basalforcings.linear.
        """
        
        if np.all(np.isnan(self.groundedice_melting_rate)):
            self.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices))
            warnings.warn('pyissm.parm.basalforcings.linear: no basalforcings.groundedice_melting_rate specified -- values set as zero.')

        return self
    
    # Marshall method for saving the basalforcings.linear parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [basalforcings.linear] parameters to a binary file.

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

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.basalforcings.model', data = 2, format = 'Integer')

        ## Write DoubleMat fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'perturbation_melting_rate', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'groundedice_melting_rate', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'geothermalflux', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofelements + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'deepwater_melting_rate', format = 'DoubleMat', mattype = 3, timeserieslength = 2, scale = 1. / md.constants.yts, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'upperwater_melting_rate', format = 'DoubleMat', mattype = 3, timeserieslength = 2, scale = 1. / md.constants.yts, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'deepwater_elevation', format = 'DoubleMat', mattype = 3, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'upperwater_elevation', format = 'DoubleMat', mattype = 3, yts = md.constants.yts)

## ------------------------------------------------------
## basalforcings.lineararma
## ------------------------------------------------------
@class_registry.register_class
class lineararma(class_registry.manage_state):
    """
    Linear ARMA basal forcings parameters class for ISSM.

    This class encapsulates the parameters for linear ARMA basal forcings in the ISSM (Ice Sheet System Model) framework.
    It defines the structure for piecewise polynomial parameters, autoregressive and moving-average coefficients, and various melting rates and elevations for different basins.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    num_basins : int, default=0
        Number of different basins [unitless].
    num_params : int, default=0
        Number of different parameters in the piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.).
    num_breaks : int, default=0
        Number of different breakpoints in the piecewise-polynomial (separating num_breaks+1 periods).
    polynomialparams : ndarray, default=np.nan
        Coefficients for the polynomial (const, trend, quadratic, etc.), dim1 for basins, dim2 for periods, dim3 for orders.
    datebreaks : ndarray, default=np.nan
        Dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr].
    ar_order : float, default=0.
        Order of the autoregressive model [unitless].
    ma_order : float, default=0.
        Order of the moving-average model [unitless].
    arma_timestep : int, default=0
        Time resolution of the ARMA model [yr].
    arlag_coefs : ndarray, default=np.nan
        Basin-specific vectors of AR lag coefficients [unitless].
    malag_coefs : ndarray, default=np.nan
        Basin-specific vectors of MA lag coefficients [unitless].
    basin_id : ndarray, default=np.nan
        Basin number assigned to each element [unitless].
    groundedice_melting_rate : ndarray, default=np.nan
        Node-specific basal melting rate for grounded ice (positive if melting) [m/yr].
    deepwater_elevation : ndarray, default=np.nan
        Basin-specific elevation of ocean deepwater [m].
    upperwater_melting_rate : ndarray, default=np.nan
        Basin-specific basal melting rate (positive if melting applied for floating ice with base >= upperwater_elevation) [m/yr].
    upperwater_elevation : ndarray, default=np.nan
        Basin-specific elevation of ocean upperwater [m].
    geothermalflux : ndarray, default=np.nan
        Node-specific geothermal heat flux [W/m^2].

    Methods
    -------
    __init__(self, other=None)
        Initializes the linear ARMA basal forcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the linear ARMA basal forcings parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.basalforcings = pyissm.model.classes.basalforcings.lineararma()
    """

    # Initialise with default parameters
    def __init__(self, other=None):
        self.num_basins = 0
        self.num_params = 0
        self.num_breaks = 0
        self.polynomialparams = np.nan
        self.datebreaks = np.nan
        self.ar_order = 0.
        self.ma_order = 0.
        self.arma_timestep = 0
        self.arlag_coefs = np.nan
        self.malag_coefs = np.nan
        self.basin_id = np.nan
        self.groundedice_melting_rate = np.nan
        self.deepwater_elevation = np.nan
        self.upperwater_melting_rate = np.nan
        self.upperwater_elevation = np.nan
        self.geothermalflux = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   basal forcings parameters:\n'

        s += '   autoregressive model is applied for deepwater_melting_rate\n'
        s += '{}\n'.format(class_utils.fielddisplay(self, 'num_basins', 'number of different basins [unitless]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'basin_id', 'basin number assigned to each element [unitless]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'num_breaks', 'number of different breakpoints in the piecewise-polynomial (separating num_breaks+1 periods)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'num_params', 'number of different parameters in the piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'polynomialparams', 'coefficients for the polynomial (const,trend,quadratic,etc.),dim1 for basins,dim2 for periods,dim3 for orders, ex: polyparams=cat(num_params,intercepts,trendlinearcoefs,trendquadraticcoefs)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'datebreaks', 'dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'ar_order', 'order of the autoregressive model [unitless]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'ma_order', 'order of the moving-average model [unitless]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'arma_timestep', 'time resolution of the ARMA model [yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'arlag_coefs', 'basin-specific vectors of AR lag coefficients [unitless]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'malag_coefs', 'basin-specific vectors of MA lag coefficients [unitless]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'deepwater_elevation', 'basin-specific elevation of ocean deepwater [m]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'upperwater_melting_rate', 'basin-specic basal melting rate (positive if melting applied for floating ice whith base >= upperwater_elevation) [m/yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'upperwater_elevation', 'basin-specific elevation of ocean upperwater [m]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'groundedice_melting_rate','node-specific basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'geothermalflux','node-specific geothermal heat flux [W/m^2]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - basalforcings.lineararma Class'

    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        if 'MasstransportAnalysis' in analyses:
            nbas = md.basalforcings.num_basins
            nprm = md.basalforcings.num_params
            nbrk = md.basalforcings.num_breaks

            class_utils.check_field(md, fieldname = "basalforcings.num_basins", scalar = True, allow_nan = False, allow_inf = False, gt = 0)
            class_utils.check_field(md, fieldname = "basalforcings.groundedice_melting_rate", timeseries = True, allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = "basalforcings.num_params", scalar = True, allow_nan = False, allow_inf = False, gt = 0)
            class_utils.check_field(md, fieldname = "basalforcings.num_breaks", scalar = True, allow_nan = False, allow_inf = False, ge = 0)

            if len(np.shape(self.deepwater_elevation)) == 1:
                self.deepwater_elevation = np.array([self.deepwater_elevation])
                self.upperwater_elevation = np.array([self.upperwater_elevation])
                self.upperwater_melting_rate = np.array([self.upperwater_melting_rate])
            if len(np.shape(self.polynomialparams)) == 1:
                self.polynomialparams = np.array([[self.polynomialparams]])

            if nbas > 1 and nbrk >= 1 and nprm > 1:
                class_utils.check_field(md, fieldname = "basalforcings.polynomialparams", allow_nan = False, allow_inf = False, size = (nbas, nbrk + 1, nprm), numel = nbas * (nbrk + 1) * nprm)
            elif nbas == 1:
                class_utils.check_field(md, fieldname = "basalforcings.polynomialparams", allow_nan = False, allow_inf = False, size = (nprm, nbrk + 1), numel = nbas * (nbrk + 1) * nprm)
            elif nbrk == 0:
                class_utils.check_field(md, fieldname = "basalforcings.polynomialparams", allow_nan = False, allow_inf = False, size = (nbas, nprm), numel = nbas * (nbrk + 1) * nprm)
            elif nprm == 1:
                class_utils.check_field(md, fieldname = "basalforcings.polynomialparams", allow_nan = False, allow_inf = False, size = (nbas, nbrk), numel = nbas * (nbrk + 1) * nprm)

            class_utils.check_field(md, fieldname = "basalforcings.deepwater_elevation", allow_nan = False, allow_inf = False, size = (1, md.basalforcings.num_basins), numel = md.basalforcings.num_basins)
            class_utils.check_field(md, fieldname = "basalforcings.upperwater_elevation", allow_nan = False, allow_inf = False, le = 0, size = (1, md.basalforcings.num_basins), numel = md.basalforcings.num_basins)
            class_utils.check_field(md, fieldname = "basalforcings.upperwater_melting_rate", allow_nan = False, allow_inf = False, ge = 0, size = (1, md.basalforcings.num_basins), numel = md.basalforcings.num_basins)
            class_utils.check_field(md, fieldname = "basalforcings.basin_id", allow_inf = False, ge = 0, le = md.basalforcings.num_basins, size = (md.mesh.numberofelements,))

            class_utils.check_field(md, fieldname = "basalforcings.ar_order", scalar = True, allow_nan = False, allow_inf = False, ge = 0)
            class_utils.check_field(md, fieldname = "basalforcings.ma_order", scalar = True, allow_nan = False, allow_inf = False, ge = 0)
            class_utils.check_field(md, fieldname = "basalforcings.arma_timestep", scalar = True, allow_nan = False, allow_inf = False, ge = md.timestepping.time_step)
            class_utils.check_field(md, fieldname = "basalforcings.arlag_coefs", allow_nan = False, allow_inf = False, size = (md.basalforcings.num_basins, md.basalforcings.ar_order))
            class_utils.check_field(md, fieldname = "basalforcings.malag_coefs", allow_nan = False, allow_inf = False, size = (md.basalforcings.num_basins, md.basalforcings.ma_order))

            if nbrk > 0:
                class_utils.check_field(md, fieldname = "basalforcings.datebreaks", allow_nan = False, allow_inf = False, size = (nbas, nbrk))
            elif np.size(md.basalforcings.datebreaks) == 0 or np.all(np.isnan(md.basalforcings.datebreaks)):
                pass
            else:
                raise RuntimeError("md.basalforcings.num_breaks is 0 but md.basalforcings.datebreaks is not empty")

        if 'BalancethicknessAnalysis' in analyses:
            raise Exception("pyissm.basalforcings.lineararma.check_consistency:: BalancethicknessAnalysis not implemented yet!")
        if 'ThermalAnalysis' in analyses and solution != 'TransientSolution' and not md.transient.isthermal:
            raise Exception("pyissm.basalforcings.lineararma.check_consistency:: ThermalAnalysis not implemented yet!")

        return md
    
    # Initialise empty fields of correct dimensions
    def initialise(self, md):
        """
        Initialise empty fields in basalforcings.lineararma.
        """

        if np.all(np.isnan(self.groundedice_melting_rate)):
            self.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices))
            warnings.warn('pyissm.model.classes.basalforcings.lineararma: no basalforcings.groundedice_melting_rate specified -- values set as 0.')

        if np.all(np.isnan(self.trend)):
            self.trend = np.zeros((1, self.num_basins)) # No trend in SMB
            warnings.warn('pyissm.model.classes.basalforcings.lineararma: no basalforcings.trend specified -- values set as 0.')

        if self.ar_order == 0:
            self.ar_order = 1 # Dummy 1 value for autoregression
            self.arlag_coefs = np.zeros((self.num_basins, self.ar_order)) # Autoregression coefficients all set to 0
            warnings.warn('pyissm.model.classes.basalforcings.lineararma: no basalforcings.ar_order specified -- order of autoregressive model set to 0.')

        if self.arma_timestep == 0:
            self.arma_timestep = md.timestepping.time_step # ARMA model has no prescribed time step
            warnings.warn('pyissm.model.classes.basalforcings.lineararma: no basalforcings.arma_timestep specified -- set to md.timestepping.time_step.')

        if np.all(np.isnan(self.arlag_coefs)):
            self.arlag_coefs = np.zeros((self.num_basins, self.ar_order)) # Autoregression model of order 0
            warnings.warn('pyissm.model.classes.basalforcings.lineararma: no basalforcings.arlag_coefs specified -- order of autoregressive model set to 0.')

        if np.all(np.isnan(self.malag_coefs)):
            self.malag_coefs = np.zeros((self.num_basins, self.ma_order)) # Moving-average model of order 0
            warnings.warn('pyissm.model.classes.basalforcings.lineararma: no basalforcings.malag_coefs specified -- order of moving-average model set to 0.')

        return self

    # Marshall method for saving the basalforcings.lineararma parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [basalforcings.lineararma] parameters to a binary file.

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

        ## Scale parameters
        ## NOTE: Scaling logic here taken from $ISSM_DIR/src/m/classes/linearbasalforcingsarma.py
        polyParams_scaled = np.copy(self.polynomialparams)
        polyParams_scaled_2d = np.zeros((self.num_basins, self.num_breaks + self.num_params))

        if(self.num_params > 1):
            # Case 3D
            if(self.num_basins > 1 and self.num_breaks + 1 > 1):
                for ii in range(self.num_params):
                    polyParams_scaled[:, :, ii] = polyParams_scaled[:, :, ii] * (1. / md.constants.yts) ** (ii + 1)
                ## Fit in 2D array
                for ii in range(self.num_params):
                    polyParams_scaled_2d[:, ii * (self.num_breaks + 1):(ii + 1) * (self.num_breaks + 1)] = 1 * polyParams_scaled[:, :, ii]
            # Case 2D and higher-order params at increasing row index #
            elif(self.num_basins == 1):
                for ii in range(self.num_params):
                    polyParams_scaled[ii, :] = polyParams_scaled[ii,:] * (1. / md.constants.yts) ** (ii + 1)
                ## Fit in row array
                for ii in range(self.num_params):
                    polyParams_scaled_2d[0, ii * (self.num_breaks + 1):(ii + 1) * (self.num_breaks + 1)] = 1 * polyParams_scaled[ii, :]
            # Case 2D and higher-order params at incrasing column index #
            elif(self.num_breaks + 1 == 1):
                for ii in range(self.num_params):
                    polyParams_scaled[:, ii] = polyParams_scaled[:, ii] * (1. / md.constants.yts) ** (ii + 1)
                # 2D array is already in correct format #
                polyParams_scaled_2d = np.copy(polyParams_scaled)
        else:
            polyParams_scaled   = polyParams_scaled * (1. / md.constants.yts)
            # 2D array is already in correct format #
            polyParams_scaled_2d = np.copy(polyParams_scaled)

        if(self.num_breaks + 1 == 1):
            dbreaks = np.zeros((self.num_basins, 1))
        else:
            dbreaks = np.copy(self.datebreaks)


        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.basalforcings.model', data = 9, format = 'Integer')

        ## Write Integer fields
        fieldnames = ['num_basins', 'num_params', 'num_breaks', 'ar_order', 'ma_order']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Integer')

        ## Write DoubleMat fields
        execute.WriteData(fid, prefix, name = 'md.basalforcings.polynomialparams', data = polyParams_scaled_2d, format = 'DoubleMat')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'arlag_coefs', format = 'DoubleMat', yts =  md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'malag_coefs', format = 'DoubleMat', yts =  md.constants.yts)
        execute.WriteData(fid, prefix, name = 'md.basalforcings.datebreaks', data = dbreaks, format = 'DoubleMat', scale = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'deepwater_elevation', format = 'DoubleMat')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'upperwater_melting_rate', format = 'DoubleMat', scale = 1. / md.constants.yts, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'upperwater_elevation', format = 'DoubleMat')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'groundedice_melting_rate', format = 'DoubleMat', scale = 1. / md.constants.yts, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'geothermalflux', format = 'DoubleMat', scale = 1. / md.constants.yts, yts = md.constants.yts)

        ## Write IntMat fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'basin_id', data = self.basin_id - 1, format = 'IntMat', mattype = 2)  # 0-indexed

        ## Write Double fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'arma_timestep', format = 'Double', scale = md.constants.yts)

## ------------------------------------------------------
## basalforcings.mismip
## ------------------------------------------------------
@class_registry.register_class
class mismip(class_registry.manage_state):
    """
    MISMIP basal forcings parameters class for ISSM.

    This class encapsulates the parameters for the MISMIP basal melt parameterization in the ISSM (Ice Sheet System Model) framework.
    It defines the basal melting rate for grounded ice, a melt rate factor, a threshold thickness for saturation of basal melting,
    an upper depth above which the melt rate is zero, and geothermal heat flux.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    groundedice_melting_rate : ndarray, default=np.nan
        Basal melting rate for grounded ice (positive if melting) [m/yr].
    meltrate_factor : float, default=0.2
        Melt-rate factor [1/yr] (sign is opposite to MISMIP+ benchmark to remain consistent with ISSM convention of positive values for melting).
    threshold_thickness : float, default=75.
        Threshold thickness for saturation of basal melting [m].
    upperdepth_melt : float, default=-100.
        Depth above which melt rate is zero [m].
    geothermalflux : float, default=np.nan
        Geothermal heat flux [W/m^2].

    Methods
    -------
    __init__(self, other=None)
        Initializes the MISMIP basal forcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the MISMIP basal forcings parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.basalforcings = pyissm.model.classes.basalforcings.mismip()
    md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices,))
    md.basalforcings.meltrate_factor = 0.2
    md.basalforcings.threshold_thickness = 75.
    md.basalforcings.upperdepth_melt = -100.
    """

    # Initialise with default parameters
    def __init__(self, other=None):
        self.groundedice_melting_rate = np.nan
        self.meltrate_factor = 0.2
        self.threshold_thickness = 75.
        self.upperdepth_melt = -100.
        self.geothermalflux = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   MISMIP + basal melt parameterization\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, "groundedice_melting_rate", "basal melting rate (positive if melting) [m / yr]"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "meltrate_factor", "Melt - rate rate factor [1 / yr] (sign is opposite to MISMIP + benchmark to remain consistent with ISSM convention of positive values for melting)"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "threshold_thickness", "Threshold thickness for saturation of basal melting [m]"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "upperdepth_melt", "Depth above which melt rate is zero [m]"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "geothermalflux", "Geothermal heat flux [W / m^2]"))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - basalforcings.mismip Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        if 'MasstransportAnalysis' in analyses and solution != 'TransientSolution' and not md.transient.ismasstransport:
            class_utils.check_field(md, fieldname = "basalforcings.groundedice_melting_rate", timeseries = True, allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = "basalforcings.meltrate_factor", scalar = True, ge = 0)
            class_utils.check_field(md, fieldname = "basalforcings.threshold_thickness", scalar = True, ge = 0)
            class_utils.check_field(md, fieldname = "basalforcings.upperdepth_melt", scalar = True, le = 0)

        if 'BalancethicknessAnalysis' in analyses:
            class_utils.check_field(md, fieldname = "basalforcings.groundedice_melting_rate", allow_nan = False, allow_inf = False, size = (md.mesh.numberofvertices, ))
            class_utils.check_field(md, fieldname = "basalforcings.meltrate_factor", scalar = True, ge = 0)
            class_utils.check_field(md, fieldname = "basalforcings.threshold_thickness", scalar = True, ge = 0)
            class_utils.check_field(md, fieldname = "basalforcings.upperdepth_melt", scalar = True, le = 0)

        if 'ThermalAnalysis' in analyses and not (solution == 'TransientSolution' and not md.transient.isthermal):
            class_utils.check_field(md, fieldname = "basalforcings.groundedice_melting_rate", timeseries = True, allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = "basalforcings.meltrate_factor", scalar = True, ge = 0)
            class_utils.check_field(md, fieldname = "basalforcings.threshold_thickness", scalar = True, ge = 0)
            class_utils.check_field(md, fieldname = "basalforcings.upperdepth_melt", scalar = True, le = 0)
            class_utils.check_field(md, fieldname = "basalforcings.geothermalflux", timeseries = True, allow_nan = False, allow_inf = False, ge = 0)

        return md
    
    # Initialise empty fields of correct dimensions
    def initialise(self, md):
        """
        Initialise empty fields in basalforcings.mismip.
        """

        if np.all(np.isnan(self.groundedice_melting_rate)):
            self.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices))
            warnings.warn('pyissm.model.classes.basalforcings.mismip: no basalforcings.groundedice_melting_rate specified -- values set as 0.')
        
        if np.all(np.isnan(self.geothermalflux)):
            self.geothermalflux = np.zeros((md.mesh.numberofvertices))
            warnings.warn('pyissm.model.classes.basalforcings.mismip: no basalforcings.geothermalflux specified -- values set as 0.')

        return self

    # Marshall method for saving the basalforcings.mismip parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [basalforcings.mismip] parameters to a binary file.

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

        ## Write warning if yts does not match this adjusted value
        if md.constants.yts != 365.2422 * 24. * 3600.:
            print('WARNING: value of yts for MISMIP + runs different from ISSM default!')

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.basalforcings.model', data = 3, format = 'Integer')

        ## Write DoubleMat fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'groundedice_melting_rate', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'geothermalflux', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        
        ## Write Double fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'meltrate_factor', format = 'Double', scale = 1. / md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'threshold_thickness', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'upperdepth_melt', format = 'Double')

## ------------------------------------------------------
## basalforcings.plume
## ------------------------------------------------------
@class_registry.register_class
class plume(class_registry.manage_state):
    """
    Plume basal forcings parameters class for ISSM.

    This class encapsulates the parameters for plume basal forcings in the ISSM (Ice Sheet System Model) framework.
    It defines the structure of the mantle plume, including its radius, depth, and position, as well as parameters related to geothermal heat flux and melting rates.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    floatingice_melting_rate : ndarray, default=np.nan
        Basal melting rate for floating ice (positive if melting) [m/yr].
    groundedice_melting_rate : ndarray, default=np.nan
        Basal melting rate for grounded ice (positive if melting) [m/yr].
    mantleconductivity : float, default=2.2
        Mantle heat conductivity [W/m^3].
    nusselt : float, default=300
        Nusselt number, ratio of mantle to plume [1].
    dtbg : float, default=0.011
        Background temperature gradient [degree/m].
    plumeradius : float, default=100000
        Radius of the mantle plume [m].
    topplumedepth : float, default=10000
        Depth of the mantle plume top below the crust [m].
    bottomplumedepth : float, default=1050000
        Depth of the mantle plume base below the crust [m].
    plumex : float, default=np.nan
        x coordinate of the center of the plume [m].
    plumey : float, default=np.nan
        y coordinate of the center of the plume [m].
    crustthickness : float, default=30000
        Thickness of the crust [m].
    uppercrustthickness : float, default=14000
        Thickness of the upper crust [m].
    uppercrustheat : float, default=1.7e-6
        Volumic heat of the upper crust [W/m^3].
    lowercrustheat : float, default=0.4e-6
        Volumic heat of the lower crust [W/m^3].

    Methods
    -------
    __init__(self, other=None)
        Initializes the plume basal forcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the plume basal forcings parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.basalforcings = pyissm.model.classes.basalforcings.plume()
    md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices,))
    md.basalforcings.floatingice_melting_rate = np.ones((md.mesh.numberofvertices,)) * 2
    """

    # Initialise with default parameters
    def __init__(self, other=None):
        self.floatingice_melting_rate = np.nan
        self.groundedice_melting_rate = np.nan
        self.mantleconductivity = 2.2
        self.nusselt = 300
        self.dtbg = 11 / 1000.
        self.plumeradius = 100000
        self.topplumedepth = 10000
        self.bottomplumedepth = 1050000
        self.plumex = np.nan
        self.plumey = np.nan
        self.crustthickness = 30000
        self.uppercrustthickness = 14000
        self.uppercrustheat = 1.7 * pow(10, -6)
        self.lowercrustheat = 0.4 * pow(10, -6)

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   mantle plume basal melt parameterization:\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'groundedice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'floatingice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'mantleconductivity', 'mantle heat conductivity [W/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'nusselt', 'nusselt number, ratio of mantle to plume [1]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'dtbg', 'background temperature gradient [degree/m]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'plumeradius', 'radius of the mantle plume [m]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'topplumedepth', 'depth of the mantle plume top below the crust [m]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'bottomplumedepth', 'depth of the mantle plume base below the crust [m]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'plumex', 'x coordinate of the center of the plume [m]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'plumey', 'y coordinate of the center of the plume [m]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'crustthickness', 'thickness of the crust [m]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'uppercrustthickness', 'thickness of the upper crust [m]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'uppercrustheat', 'volumic heat of the upper crust [w/m^3]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'lowercrustheat', 'volumic heat of the lowercrust [w/m^3]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - basalforcings.plume Class'
        return s
    
    # Check model consistency
    def checkconsistency(self, md, solution, analyses):
        if 'MasstransportAnalysis' in analyses and not (solution == 'TransientSolution' and md.transient.ismasstransport == 0):
            class_utils.check_field(md, fieldname = "basalforcings.groundedice_melting_rate", timeseries = True, allow_nan = False)
            class_utils.check_field(md, fieldname = "basalforcings.floatingice_melting_rate", timeseries = True, allow_nan = False)

        if 'BalancethicknessAnalysis' in analyses:
            class_utils.check_field(md, fieldname = "basalforcings.groundedice_melting_rate", size = (md.mesh.numberofvertices, ), allow_nan = False)
            class_utils.check_field(md, fieldname = "basalforcings.floatingice_melting_rate", size = (md.mesh.numberofvertices, ), allow_nan = False)

        if 'ThermalAnalysis' in analyses and not (solution == 'TransientSolution' and md.transient.isthermal == 0):
            class_utils.check_field(md, fieldname = "basalforcings.groundedice_melting_rate", timeseries = True, allow_nan = False)
            class_utils.check_field(md, fieldname = "basalforcings.floatingice_melting_rate", timeseries = True, allow_nan = False)
            class_utils.check_field(md, fieldname = "basalforcings.mantleconductivity", scalar = True, ge = 0)
            class_utils.check_field(md, fieldname = "basalforcings.nusselt", scalar = True, gt = 0)
            class_utils.check_field(md, fieldname = "basalforcings.dtbg", scalar = True, gt = 0)
            class_utils.check_field(md, fieldname = "basalforcings.topplumedepth", scalar = True, gt = 0)
            class_utils.check_field(md, fieldname = "basalforcings.bottomplumedepth", scalar = True, gt = 0)
            class_utils.check_field(md, fieldname = "basalforcings.plumex", scalar = True)
            class_utils.check_field(md, fieldname = "basalforcings.plumey", scalar = True)
            class_utils.check_field(md, fieldname = "basalforcings.crustthickness", scalar = True, gt = 0)
            class_utils.check_field(md, fieldname = "basalforcings.uppercrustthickness", scalar = True, gt = 0)
            class_utils.check_field(md, fieldname = "basalforcings.uppercrustheat", scalar = True, gt = 0)
            class_utils.check_field(md, fieldname = "basalforcings.lowercrustheat", scalar = True, gt = 0)

        return md

    # Initialise empty fields of correct dimensions
    def initialise(self, md):
        """
        Initialise empty fields in basalforcings.plume.
        """

        if np.all(np.isnan(self.groundedice_melting_rate)):
            self.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices,))
            warnings.warn('pyissm.model.classes.basalforcings.plume: no groundedice_melting_rate specified -- values set as 0')

        if np.all(np.isnan(self.floatingice_melting_rate)):
            self.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices,))
            warnings.warn('pyissm.model.classes.basalforcings.plume: no floatingice_melting_rate specified -- values set as 0')

        return self

    # Marshall method for saving the basalforcings.plume parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [basalforcings.plume] parameters to a binary file.

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

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.basalforcings.model', data = 4, format = 'Integer')

        ## Write DoubleMat fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'groundedice_melting_rate', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'floatingice_melting_rate', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        
        ## Write Double fields
        fieldnames = ['mantleconductivity', 'nusselt', 'dtbg', 'plumeradius', 'topplumedepth',
                      'bottomplumedepth', 'plumex', 'plumey', 'crustthickness', 'uppercrustthickness',
                      'uppercrustheat', 'lowercrustheat']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Double')

## ------------------------------------------------------
## basalforcings.spatiallinear
## ------------------------------------------------------
@class_registry.register_class
class spatiallinear(class_registry.manage_state):
    """
    Spatial linear basal forcings parameters class for ISSM.

    This class encapsulates the parameters for spatial linear basal forcings in the ISSM (Ice Sheet System Model) framework.
    It defines the melting rates for grounded ice, deepwater, and upperwater, as well as geothermal heat flux and perturbation melting rate.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    groundedice_melting_rate : ndarray, default=np.nan
        Basal melting rate for grounded ice (positive if melting) [m/yr].
    deepwater_melting_rate : ndarray, default=np.nan
        Basal melting rate applied for floating ice with base < deepwater_elevation [m/yr].
    deepwater_elevation : ndarray, default=np.nan
        Elevation threshold for deepwater melting rate [m].
    upperwater_melting_rate : ndarray, default=np.nan
        Basal melting rate applied for floating ice with base >= upperwater_elevation [m/yr].
    upperwater_elevation : ndarray, default=np.nan
        Elevation threshold for upperwater melting rate [m].
    geothermalflux : ndarray, default=np.nan
        Geothermal heat flux [W/m^2].
    perturbation_melting_rate : ndarray, default=np.nan
        Basal melting rate perturbation added to computed melting rate (positive if melting) [m/yr].

    Methods
    -------
    __init__(self, other=None)
        Initializes the spatial linear basal forcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the spatial linear basal forcings parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.basalforcings = pyissm.model.classes.basalforcings.spatiallinear()
    """

    # Initialise with default parameters
    def __init__(self, other=None):
        self.groundedice_melting_rate = np.nan
        self.deepwater_melting_rate = np.nan
        self.deepwater_elevation = np.nan
        self.upperwater_melting_rate = np.nan
        self.upperwater_elevation = np.nan
        self.geothermalflux = np.nan
        self.perturbation_melting_rate = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   spatial linear basal forcings parameters:\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'groundedice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'deepwater_melting_rate', 'basal melting rate (positive if melting applied for floating ice whith base < deepwater_elevation) [m/yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'deepwater_elevation', 'elevation of ocean deepwater [m]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'upperwater_melting_rate', 'basal melting rate (positive if melting applied for floating ice whith base >= upperwater_elevation) [m/yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'upperwater_elevation', 'elevation of ocean upperwater [m]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'perturbation_melting_rate', 'basal melting rate perturbation added to computed melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'geothermalflux', 'geothermal heat flux [W/m^2]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - basalforcings.spatiallinear Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        if not np.all(np.isnan(self.perturbation_melting_rate)):
            class_utils.check_field(md, fieldname = "basalforcings.perturbation_melting_rate", timeseries = True, allow_nan = False)

        if 'MasstransportAnalysis' in analyses and not solution == 'TransientSolution' and not md.transient.ismasstransport:
            class_utils.check_field(md, fieldname = "basalforcings.groundedice_melting_rate", timeseries = True, allow_nan = False)
            class_utils.check_field(md, fieldname = "basalforcings.deepwater_melting_rate", timeseries = True, ge = 0, allow_nan = False)
            class_utils.check_field(md, fieldname = "basalforcings.deepwater_elevation", timeseries = True, allow_nan = False)
            class_utils.check_field(md, fieldname = "basalforcings.upperwater_melting_rate", timeseries = True, ge = 0, allow_nan = False)
            class_utils.check_field(md, fieldname = "basalforcings.upperwater_elevation", timeseries = True, lt = 0, allow_nan = False)

        if 'BalancethicknessAnalysis' in analyses:
            raise Exception("pyissm.model.classes.basalforcings.spatiallinear.check_consistency:: BalancethicknessAnalysis not implemented yet!")

        if 'ThermalAnalysis' in analyses and not solution == 'TransientSolution' and not md.transient.isthermal:
            raise Exception("pyissm.model.classes.basalforcings.spatiallinear.check_consistency:: ThermalAnalysis not implemented yet!")

        return md
        
    # Initialise empty fields of correct dimensions
    def initialise(self, md):
        """
        Initialise empty fields in basalforcings.spatiallinear.
        """

        if np.all(np.isnan(self.groundedice_melting_rate)):
            self.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices,))
            warnings.warn('pyissm.model.classes.basalforcings.spatiallinear: no groundedice_melting_rate specified -- values set as 0')

        return self

    # Marshall method for saving the basalforcings.spatiallinear parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [basalforcings.spatiallinear] parameters to a binary file.

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

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.basalforcings.model', data = 6, format = 'Integer')

        ## Write DoubleMat fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'groundedice_melting_rate', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'geothermalflux', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'deepwater_melting_rate', format = 'DoubleMat', scale = 1. / md.constants.yts, mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'deepwater_elevation', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'upperwater_melting_rate', format = 'DoubleMat', scale = 1. / md.constants.yts, mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'upperwater_elevation', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'perturbation_melting_rate', format = 'DoubleMat',  scale = 1. / md.constants.yts, mattype = 1)
