"""
Frontal Forcings classes for ISSM.
"""

import numpy as np
import warnings
from pyissm.model.classes import class_utils, class_registry
from pyissm.model import execute, mesh

## ------------------------------------------------------
## frontalforcings.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    """
    Default frontalforcings class for ISSM.

    This class contains the default parameters for frontal forcings in the ISSM framework.
    It defines the main frontal forcing-related parameters.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    meltingrate : :class:`numpy.ndarray`, default=np.nan
        Melting rate at given location [m/a].
    ablationrate : :class:`numpy.ndarray`, default=np.nan
        Frontal ablation rate at given location [m/a], it contains both calving and melting.

    Examples
    --------
    .. code-block:: python
    
        >>> md.frontalforcings = pyissm.model.classes.frontalforcings.default()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.meltingrate = np.nan
        self.ablationrate = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Frontalforcings parameters:\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'meltingrate', 'melting rate at given location [m/a]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'ablationrate', 'frontal ablation rate at given location [m/a], it contains both calving and melting'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - frontalforcings.default Class'
        return s
    
    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude frontalforcings.default fields to 3D
        """
        self.meltingrate = mesh.project_3d(md, vector = self.meltingrate, type = 'node')
        self.ablationrate = mesh.project_3d(md, vector = self.ablationrate, type = 'node')
            
        return self

    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [frontalforcings.default] parameters.

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

        # Early return if not transient movingfront analysis
        if (solution != 'TransientSolution') or (not md.transient.ismovingfront):
            return md

        class_utils.check_field(md, fieldname = 'frontalforcings.meltingrate', timeseries = True, ge = 0, allow_nan = False, allow_inf = False)
        if not np.isnan(md.frontalforcings.ablationrate):
            class_utils.check_field(md, fieldname = 'frontalforcings.ablationrate', timeseries = True, allow_nan = False, allow_inf = False)
            
        return md
    
    # Marshall method for saving the frontalforcings.default parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [frontalforcings.default] parameters to a binary file.

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
        execute.WriteData(fid, prefix, name = 'md.frontalforcings.parameterization', data = 1, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'meltingrate', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts, scale = 1. / md.constants.yts)

        ## Write conditional field
        if not np.isnan(self.ablationrate).all():
            execute.WriteData(fid, prefix, obj = self, fieldname = 'ablationrate', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts, scale = 1. / md.constants.yts)

## ------------------------------------------------------
## frontalforcings.rignot
## ------------------------------------------------------
@class_registry.register_class
class rignot(class_registry.manage_state):
    """
    Rignot frontalforcings class for ISSM.

    This class contains the parameters for frontal forcings based on the Rignot methodology in the ISSM framework.
    It defines the main frontal forcing-related parameters specific to the Rignot approach.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    basin_id : :class:`numpy.ndarray`, default=np.nan
        Basin ID for elements.
    num_basins : :class:`int`, default=0
        Number of basins.
    subglacial_discharge : :class:`numpy.ndarray`, default=np.nan
        Sum of subglacial discharge for each basin [m/d].
    thermalforcing : :class:`numpy.ndarray`, default=np.nan
        Thermal forcing [°C].

    Examples
    --------
    .. code-block:: python
    
        >>> md.frontalforcings = pyissm.model.classes.frontalforcings.rignot()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.basin_id = np.nan
        self.num_basins = 0
        self.subglacial_discharge = np.nan
        self.thermalforcing = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Frontalforcings parameters:\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'basin_id', 'basin ID for elements'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'num_basins', 'number of basins'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'subglacial_discharge', 'sum of subglacial discharge for each basin [m/d]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'thermalforcing', 'thermal forcing [∘C]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - frontalforcings.rignot Class'
        return s
    
    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude frontalforcings.rignot fields to 3D
        """
        warnings.warn('pyissm.model.classes.frontalforcings.rignot.extrude: 3D extrusion not implemented for frontalforcings.rignot. Returning unchanged (2D) frontalforcing fields.')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [frontalforcings.rignot] parameters.

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

        # Early return
        if (solution != 'TransientSolution') or (not md.transient.ismovingfront):
            return md

        class_utils.check_field(md, fieldname = 'frontalforcings.num_basins', scalar = True, gt = 0, allow_nan = True, allow_inf = True)
        class_utils.check_field(md, fieldname = 'frontalforcings.basin_id', ge = 0, le = md.frontalforcings.num_basins, size = (md.mesh.numberofelements, ), allow_inf = False)
        class_utils.check_field(md, fieldname = 'frontalforcings.subglacial_discharge', timeseries = True, ge = 0, allow_nan = False, allow_inf = True)
        class_utils.check_field(md, fieldname = 'frontalforcings.thermalforcing', timeseries = True, allow_nan = False, allow_inf = False)

        return md
    
    # Marshall method for saving the frontalforcings.rignot parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [frontalforcings.rignot] parameters to a binary file.

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
        execute.WriteData(fid, prefix, name = 'md.frontalforcings.parameterization', data = 2, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, name = 'md.frontalforcings.basin_id', data = self.basin_id - 1, format = 'IntMat', mattype = 2) # 0-indexed
        execute.WriteData(fid, prefix, obj = self, fieldname = 'num_basins', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'subglacial_discharge', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'thermalforcing', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

## ------------------------------------------------------
## frontalforcings.rignotarma
## ------------------------------------------------------
@class_registry.register_class
class rignotarma(class_registry.manage_state):
    """
    RignotARMA frontalforcings class for ISSM.

    This class contains the parameters for frontal forcings based on the Rignot 
    methodology with ARMA (AutoRegressive Moving Average) modeling in the ISSM framework.
    It defines the main frontal forcing-related parameters specific to the RignotARMA approach,
    including polynomial trends, breakpoints, ARMA coefficients, and subglacial discharge modeling.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    num_basins : :class:`int`, default=0
        Number of different basins.
    num_params : :class:`int`, default=0
        Number of different parameters in the piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.).
    num_breaks : :class:`int`, default=0
        Number of different breakpoints in the piecewise-polynomial (separating num_breaks+1 periods).
    polynomialparams : :class:`numpy.ndarray`, default=np.nan
        Coefficients for the polynomial (const, trend, quadratic, etc.), dim1 for basins, dim2 for periods, dim3 for orders.
    datebreaks : :class:`numpy.ndarray`, default=np.nan
        Dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr].
    ar_order : :class:`int`, default=0
        Order of the autoregressive model.
    ma_order : :class:`int`, default=0
        Order of the moving-average model.
    arma_timestep : :class:`int`, default=0
        Time resolution of the ARMA model [yr].
    arlag_coefs : :class:`numpy.ndarray`, default=np.nan
        Basin-specific vectors of AR lag coefficients.
    malag_coefs : :class:`numpy.ndarray`, default=np.nan
        Basin-specific vectors of MA lag coefficients.
    monthlyvals_intercepts : :class:`numpy.ndarray`, default=np.nan
        Monthly intercept values for each basin.
    monthlyvals_trends : :class:`numpy.ndarray`, default=np.nan
        Monthly trend values for each basin.
    monthlyvals_numbreaks : :class:`int`, default=0
        Number of breakpoints for monthly values.
    monthlyvals_datebreaks : :class:`numpy.ndarray`, default=np.nan
        Dates at which the monthly value breakpoints occur.
    basin_id : :class:`numpy.ndarray`, default=np.nan
        Basin number assigned to each element.
    subglacial_discharge : :class:`numpy.ndarray`, default=np.nan
        Sum of subglacial discharge for each basin [m/d].
    isdischargearma : :class:`int`, default=0
        Whether an ARMA model is also used for the subglacial discharge (if 0: subglacial_discharge is used, if 1: sd_ parameters are used).
    sd_ar_order : :class:`int`, default=0
        Order of the subglacial discharge autoregressive model.
    sd_ma_order : :class:`int`, default=0
        Order of the subglacial discharge moving-average model.
    sd_arma_timestep : :class:`int`, default=0
        Time resolution of the subglacial discharge ARMA model [yr].
    sd_arlag_coefs : :class:`numpy.ndarray`, default=np.nan
        Basin-specific vectors of AR lag coefficients for subglacial discharge.
    sd_malag_coefs : :class:`numpy.ndarray`, default=np.nan
        Basin-specific vectors of MA lag coefficients for subglacial discharge.
    sd_monthlyfrac : :class:`numpy.ndarray`, default=np.nan
        Basin-specific vectors of 12 values with fraction of the annual discharge occurring every month.
    sd_num_breaks : :class:`int`, default=0
        Number of different breakpoints in the subglacial discharge piecewise-polynomial (separating sd_num_breaks+1 periods).
    sd_num_params : :class:`int`, default=0
        Number of different parameters in the subglacial discharge piecewise-polynomial.
    sd_polynomialparams : :class:`numpy.ndarray`, default=np.nan
        Coefficients for the subglacial discharge polynomial (const, trend, quadratic, etc.).
    sd_datebreaks : :class:`numpy.ndarray`, default=np.nan
        Dates at which the breakpoints in the subglacial discharge piecewise polynomial occur (1 row per basin) [yr].

    Examples
    --------
    .. code-block:: python
    
        >>> md.frontalforcings = pyissm.model.classes.frontalforcings.rignotarma()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.num_basins = 0
        self.num_params = 0
        self.num_breaks = 0
        self.polynomialparams = np.nan
        self.datebreaks       = np.nan
        self.ar_order = 0
        self.ma_order = 0
        self.arma_timestep = 0
        self.arlag_coefs = np.nan
        self.malag_coefs = np.nan
        self.monthlyvals_intercepts = np.nan
        self.monthlyvals_trends = np.nan
        self.monthlyvals_numbreaks = 0
        self.monthlyvals_datebreaks = np.nan
        self.basin_id = np.nan
        self.subglacial_discharge = np.nan
        self.isdischargearma = 0
        self.sd_ar_order = 0.
        self.sd_ma_order = 0.
        self.sd_arma_timestep = 0
        self.sd_arlag_coefs = np.nan
        self.sd_malag_coefs = np.nan
        self.sd_monthlyfrac = np.nan
        self.sd_num_breaks  = 0
        self.sd_num_params  = 0
        self.sd_polynomialparams = np.nan
        self.sd_datebreaks = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Frontalforcings parameters:\n'

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
        s += '{}\n'.format(class_utils.fielddisplay(self, 'isdischargearma','whether an ARMA model is also used for the subglacial discharge (if 0: subglacial_discharge is used, if 1: sd_ parameters are used)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'subglacial_discharge', 'sum of subglacial discharge for each basin [m/d]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'sd_ar_order','order of the subglacial discharge autoregressive model [unitless]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'sd_ma_order','order of the subglacial discharge moving-average model [unitless]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'sd_arma_timestep','time resolution of the subglacial discharge autoregressive model [yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'sd_arlag_coefs','basin-specific vectors of AR lag coefficients for subglacial discharge [unitless]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'sd_malag_coefs','basin-specific vectors of MA lag coefficients for subglacial discharge [unitless]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'sd_monthlyfrac','basin-specific vectors of 12 values with fraction of the annual discharge occuring every month [unitless]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'sd_num_params','number of different parameters in the subglacial discharge piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'sd_num_breaks','number of different breakpoints in the subglacial discharge piecewise-polynomial (separating sd_num_breaks+1 periods)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'sd_datebreaks','dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'sd_polynomialparams','coefficients for the sd_polynomial (const,trend,quadratic,etc.),dim1 for basins,dim2 for periods,dim3 for orders'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - frontalforcings.rignotarma Class'
        return s
    
    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude frontalforcings.rignotarma fields to 3D
        """
        warnings.warn('pyissm.model.classes.frontalforcings.rignotarma.extrude: 3D extrusion not implemented for frontalforcings.rignotarma. Returning unchanged (2D) frontalforcing fields.')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [frontalforcings.rignotarma] parameters.

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

        ## NOTE: Logic checks here taken from $ISSM_DIR/src/m/classes/frontalforcingsrignotarma.py

        # Early return if not transient movingfront analysis
        if not (solution == 'TransientSolution') or not md.transient.ismovingfront:
            return md

        nbas  = md.frontalforcings.num_basins
        nprm  = md.frontalforcings.num_params
        nbrk  = md.frontalforcings.num_breaks
        nMbrk = md.frontalforcings.monthlyvals_numbreaks

        class_utils.check_field(md, fieldname = 'frontalforcings.num_basins', scalar = True, gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'frontalforcings.num_params', scalar = True, gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'frontalforcings.num_breaks', scalar = True, ge = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'frontalforcings.basin_id', ge = 0, le = md.frontalforcings.num_basins, size = (md.mesh.numberofelements, ), allow_inf = False)
        class_utils.check_field(md, fieldname = 'frontalforcings.subglacial_discharge', timeseries = True, ge = 0, allow_nan = False, allow_inf = False)
        
        if len(np.shape(self.polynomialparams)) == 1:
            self.polynomialparams = np.array([[self.polynomialparams]])
        if(nbas>1 and nbrk>=1 and nprm>1):
            class_utils.check_field(md, fieldname = 'frontalforcings.polynomialparams', size = (nbas, nbrk+1, nprm), numel = nbas*(nbrk+1)*nprm, allow_nan = False, allow_inf = False)
        elif(nbas==1):
            class_utils.check_field(md, fieldname = 'frontalforcings.polynomialparams', size = (nprm, nbrk+1), numel = nbas*(nbrk+1)*nprm, allow_nan = False, allow_inf = False)
        elif(nbrk==0):
            class_utils.check_field(md, fieldname = 'frontalforcings.polynomialparams', size = (nbas, nprm), numel = nbas*(nbrk+1)*nprm, allow_nan = False, allow_inf = False)
        elif(nprm==1):
            class_utils.check_field(md, fieldname = 'frontalforcings.polynomialparams', size = (nbas, nbrk), numel = nbas*(nbrk+1)*nprm, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'frontalforcings.ar_order', scalar = True, ge = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'frontalforcings.ma_order', scalar = True, ge = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'frontalforcings.arma_timestep', scalar = True, ge = md.timestepping.time_step, allow_nan = False, allow_inf = False) # ARMA time step cannot be finer than ISSM timestep
        class_utils.check_field(md, fieldname = 'frontalforcings.arlag_coefs', size = (md.frontalforcings.num_basins, md.frontalforcings.ar_order), allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'frontalforcings.malag_coefs', size = (md.frontalforcings.num_basins, md.frontalforcings.ma_order), allow_nan = False, allow_inf = False)
        if(nbrk>0):
            class_utils.check_field(md, fieldname = 'frontalforcings.datebreaks', size = (nbas, nbrk), allow_nan = False, allow_inf = False)
        elif(np.size(md.frontalforcings.datebreaks)==0 or np.all(np.isnan(md.frontalforcings.datebreaks))):
            pass
        else:
            raise RuntimeError('pyissm.model.classes.frontalforcings.rignotarma.check_consistency: md.frontalforcings.num_breaks is 0 but md.frontalforcings.datebreaks is not empty')
        

        ### Check if some monthly forcings are provided ###
        if(np.all(np.isnan(md.frontalforcings.monthlyvals_intercepts))==False or np.all(np.isnan(md.frontalforcings.monthlyvals_trends))==False or np.all(np.isnan(md.frontalforcings.monthlyvals_datebreaks))==False):
            isMonthly = True
        else:
            isMonthly = False
        if(np.all(np.isnan(md.frontalforcings.monthlyvals_datebreaks))):
            isMonthlyTrend = True
        else:
            isMonthlyTrend = False
        if(isMonthly):
            class_utils.check_field(md, fieldname = 'frontalforcings.monthlyvals_numbreaks', scalar = True, ge = 0, allow_nan = False, allow_inf = False)
            if(nbas>1 and nMbrk>=1):
                class_utils.check_field(md, fieldname = 'frontalforcings.monthlyvals_intercepts', size = (nbas, 12, nMbrk+1), numel = nbas*(nMbrk+1)*12, allow_nan = False, allow_inf = False)
                if(isMonthlyTrend):
                    class_utils.check_field(md, fieldname = 'frontalforcings.monthlyvals_trends', size = (nbas, 12, nMbrk+1), numel = nbas*(nMbrk+1)*12, allow_nan = False, allow_inf = False)
            elif(nbas==1):
               class_utils.check_field(md, fieldname = 'frontalforcings.monthlyvals_intercepts', size = (nMbrk+1, 12), numel = nbas*(nMbrk+1)*12, allow_nan = False, allow_inf = False)
               if(isMonthlyTrend):
                  class_utils.check_field(md, fieldname = 'frontalforcings.monthlyvals_trends', size = (nMbrk+1, 12), numel = nbas*(nMbrk+1)*12, allow_nan = False, allow_inf = False)
            elif(nMbrk==0):
               class_utils.check_field(md, fieldname = 'frontalforcings.monthlyvals_intercepts', size = (nbas, 12), numel = nbas*(nMbrk+1)+12, allow_nan = False, allow_inf = False)
               if(isMonthlyTrend):
                  class_utils.check_field(md, fieldname = 'frontalforcings.monthlyvals_trends', size = (nbas, 12), numel = nbas*(nMbrk+1)*12, allow_nan = False, allow_inf = False)
        if(nMbrk>0):
            class_utils.check_field(md, fieldname = 'frontalforcings.monthlyvals_datebreaks', size = (nbas,nMbrk), allow_nan = False, allow_inf = False)
        elif(np.size(md.frontalforcings.monthlyvals_datebreaks)==0 or np.all(np.isnan(md.frontalforcings.monthlyvals_datebreaks))):
            pass
        else:
            raise RuntimeError('md.frontalforcings.monthlyvals_numbreaks is 0 but md.frontalforcings.monthlyvals_datebreaks is not empty')

        ### Chacking subglacial discharge ###
        class_utils.check_field(md, fieldname = 'frontalforcings.isdischargearma', scalar = True, values = [0, 1])
        if(self.isdischargearma==0):
            class_utils.check_field(md, fieldname = 'frontalforcings.subglacial_discharge', timeseries = True, ge = 0, allow_nan = False, allow_inf = False)
        else:
            sdnbrk  = md.frontalforcings.sd_num_breaks
            sdnprm  = md.frontalforcings.sd_num_params
            class_utils.check_field(md, fieldname = 'frontalforcings.sd_ar_order', scalar = True, ge = 0, allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = 'frontalforcings.sd_ma_order',scalar = True, ge = 0, allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = 'frontalforcings.sd_arma_timestep', scalar = True, ge = max(1, md.timestepping.time_step), allow_nan = False, allow_inf = False) #ARMA time step cannot be finer than ISSM timestep and annual timestep
            class_utils.check_field(md, fieldname = 'frontalforcings.sd_arlag_coefs', size = (md.frontalforcings.num_basins, md.frontalforcings.sd_ar_order), allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = 'frontalforcings.sd_malag_coefs', size = (md.frontalforcings.num_basins, md.frontalforcings.sd_ma_order), allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = 'frontalforcings.sd_monthlyfrac', size = (md.frontalforcings.num_basins, 12), allow_nan = False, allow_inf = False)
            if(np.any(abs(np.sum(self.sd_monthlyfrac,axis=1)-1)>1e-3)):
                raise RuntimeError('the 12 entries for each basin of md.frontalforcings.sd_monthlyfrac should add up to 1')
            class_utils.check_field(md, fieldname = 'frontalforcings.sd_num_params', scalar = True, gt = 0, allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = 'frontalforcings.sd_num_breaks', scalar = True, ge = 0, allow_nan = False, allow_inf = False)
            if len(np.shape(self.sd_polynomialparams)) == 1:
                self.sd_polynomialparams = np.array([[self.sd_polynomialparams]])
            if(nbas>1 and sdnbrk>=1 and sdnprm>1):
                class_utils.check_field(md, fieldname = 'frontalforcings.sd_polynomialparams', size = (nbas, sdnbrk+1, sdnprm), numel = nbas*(sdnbrk+1)*sdnprm, allow_nan = False, allow_inf = False)
            elif(nbas==1):
                class_utils.check_field(md, fieldname = 'frontalforcings.sd_polynomialparams', size = (nprm, nbrk+1), numel = nbas*(sdnbrk+1)*sdnprm, allow_nan = False, allow_inf = False)
            elif(sdnbrk==0):
                class_utils.check_field(md, fieldname = 'frontalforcings.sd_polynomialparams', size = (nbas, sdnprm), numel = nbas*(sdnbrk+1)*sdnprm, allow_nan = False, allow_inf = False)
            elif(sdnprm==1):
                class_utils.check_field(md, fieldname = 'frontalforcings.sd_polynomialparams', size = (nbas, sdnbrk), numel = nbas*(sdnbrk+1)*sdnprm, allow_nan = False, allow_inf = False)
            if(sdnbrk>0):
                class_utils.check_field(md, fieldname = 'frontalforcings.sd_datebreaks', size = (nbas, sdnbrk), allow_nan = False, allow_inf = False)
            elif(np.size(md.frontalforcings.sd_datebreaks)==0 or np.all(np.isnan(md.frontalforcings.sd_datebreaks))):
                pass
            else:
                raise RuntimeError('md.frontalforcings.sd_num_breaks is 0 but md.frontalforcings.sd_datebreaks is not empty')

        return md
    
    # Marshall method for saving the frontalforcings.rignotarma parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [frontalforcings.rignotarma] parameters to a binary file.

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

        ## Scale parameters
        ## NOTE: Scaling logic here taken from $ISSM_DIR/src/m/classes/frontalforcingsrignotarma.py
        nper = self.num_breaks + 1
        polyParams_scaled = np.copy(self.polynomialparams)
        polyParams_scaled_2d = np.zeros((self.num_basins, nper * self.num_params))
        

        if(self.num_params > 1):
            # Case 3D #
            if(self.num_basins > 1 and nper > 1):
                for ii in range(self.num_params):
                    polyParams_scaled[:,:,ii] = polyParams_scaled[:, :, ii] * (1 / md.constants.yts) ** ii
                # Fit in 2D array #
                for ii in range(self.num_params):
                    polyParams_scaled_2d[:, ii * nper:(ii + 1) * nper] = 1 * polyParams_scaled[:, :, ii]
            # Case 2D and higher-order params at increasing row index #
            elif(self.num_basins == 1):
                for ii in range(self.num_params):
                    polyParams_scaled[ii, :] = polyParams_scaled[ii, :] * (1 / md.constants.yts) ** ii
                # Fit in row array #
                for ii in range(self.num_params):
                    polyParams_scaled_2d[0, ii * nper:(ii + 1) * nper] = 1 * polyParams_scaled[ii, :]
            # Case 2D and higher-order params at increasing column index #
            elif(nper == 1):
                for ii in range(self.num_params):
                    polyParams_scaled[:, ii] = polyParams_scaled[:, ii] * (1 / md.constants.yts) ** ii
                # 2D array is already in correct format #
                polyParams_scaled_2d = np.copy(polyParams_scaled)
        else:
            # 2D array is already in correct format and no need for scaling #
            polyParams_scaled_2d = np.copy(polyParams_scaled)
        if(nper == 1):
            dbreaks = np.zeros((self.num_basins, 1))
        else:
            dbreaks = np.copy(self.datebreaks)

        ### Deal with monthly effects ###
        if(np.any(np.isnan(self.monthlyvals_intercepts))):
            interceptsM = np.zeros((self.num_basins,12)) #monthly intercepts not provided, set to 0
            trendsM     = np.zeros((self.num_basins,12)) #set monthly trends also to 0
        else:
            interceptsM3d = self.monthlyvals_intercepts
            if(np.any(np.isnan(self.monthlyvals_trends))):
                trendsM3d = 0*interceptsM3d #monthly trends not provided, set to 0
            else:
                trendsM3d = self.monthlyvals_trends
        # Create 2D arrays from 3D arrays if needed #
        if((self.monthlyvals_numbreaks + 1) > 1 and np.all(np.isnan(self.monthlyvals_intercepts))==False):
            interceptsM = np.zeros((self.num_basins, 12 * (self.monthlyvals_numbreaks + 1))) 
            trendsM     = np.zeros((self.num_basins, 12 * (self.monthlyvals_numbreaks + 1)))
            for ii in range((self.monthlyvals_numbreaks + 1)):
                interceptsM[:, ii * 12 : (ii + 1) * 12] = 1 * interceptsM3d[:,:,ii]
                trendsM[:, ii * 12 : (ii + 1) * 12] = 1 * trendsM3d[:,:,ii]
        elif((self.monthlyvals_numbreaks + 1) == 1 and np.all(np.isnan(self.monthlyvals_intercepts)) == False):
            interceptsM = 1 * interceptsM3d
            trendsM     = 1 * trendsM3d
        if((self.monthlyvals_numbreaks + 1) == 1):
            dMbreaks = np.zeros((self.num_basins, 1))
        else:
            dMbreaks = np.copy(self.monthlyvals_datebreaks)

        ### Deal with the subglacial discharge polynomial ###
        if(self.isdischargearma):
            sd_nper = self.sd_num_breaks + 1
            sdpolyParams_scaled   = np.copy(self.sd_polynomialparams)
            sdpolyParams_scaled_2d = np.zeros((self.num_basins, sd_nper * self.sd_num_params))
            if(self.sd_num_params > 1):
                # Case 3D #
                if(self.num_basins>1 and sd_nper > 1):
                    for ii in range(self.sd_num_params):
                        sdpolyParams_scaled[:, :, ii] = sdpolyParams_scaled[:, :, ii] * (1 / md.constants.yts) ** ii
                    # Fit in 2D array #
                    for ii in range(self.sd_num_params):
                        sdpolyParams_scaled_2d[:, ii * sd_nper : (ii + 1) * sd_nper] = 1 * sdpolyParams_scaled[:, :, ii]
                # Case 2D and higher-order params at increasing row index #
                elif(self.num_basins == 1):
                    for ii in range(self.sd_num_params):
                        sdpolyParams_scaled[ii, :] = sdpolyParams_scaled[ii, :] * (1 / md.constants.yts) ** ii
                    # Fit in row array #
                    for ii in range(self.num_params):
                        sdpolyParams_scaled_2d[0, ii * sd_nper : (ii + 1) * sd_nper] = 1 * sdpolyParams_scaled[ii, :]
                # Case 2D and higher-order params at incrasing column index #
                elif(sd_nper == 1):
                    for ii in range(self.sd_num_params):
                        sdpolyParams_scaled[:, ii] = sdpolyParams_scaled[:, ii] * (1 / md.constants.yts) ** ii
                    # 2D array is already in correct format #
                    sdpolyParams_scaled_2d = np.copy(sdpolyParams_scaled)
            else:
                # 2D array is already in correct format and no need for scaling #
                sdpolyParams_scaled_2d = np.copy(sdpolyParams_scaled)
            if(sd_nper == 1):
                sd_dbreaks = np.zeros((self.num_basins, 1))
            else:
                sd_dbreaks = np.copy(self.sd_datebreaks)


        ## Write headers to file
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.frontalforcings.parameterization', data = 3, format = 'Integer')

        ## Write Integer fields
        fieldnames = ['num_basins', 'num_breaks', 'num_params', 'ar_order', 'ma_order', 'monthlyvals_numbreaks']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Integer')

        ## Write DoubleMat fields
        execute.WriteData(fid, prefix, name = 'md.frontalforcings.polynomialparams', data = polyParams_scaled_2d, format = 'DoubleMat')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'arlag_coefs', format = 'DoubleMat', yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'malag_coefs', format = 'DoubleMat', yts = md.constants.yts)
        execute.WriteData(fid, prefix, name = 'md.frontalforcings.datebreaks', data = dbreaks, format = 'DoubleMat', scale = md.constants.yts)
        execute.WriteData(fid, prefix, name = 'md.frontalforcings.monthlyvals_datebreaks', data = dMbreaks, format = 'DoubleMat', scale = md.constants.yts)
        execute.WriteData(fid, prefix, name = 'md.frontalforcings.monthlyvals_intercepts', data = interceptsM, format = 'DoubleMat')
        execute.WriteData(fid, prefix, name = 'md.frontalforcings.monthlyvals_trends', data = trendsM, format = 'DoubleMat', scale = 1. / md.constants.yts)

        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'isdischargearma', format = 'Boolean')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'arma_timestep', format = 'Double', scale = md.constants.yts)
        execute.WriteData(fid, prefix, name = 'md.frontalforcings.basin_id', data =  self.basin_id - 1, format = 'IntMat', mattype = 2)  # 0-indexed

        ## Write conditional fields
        if(self.isdischargearma == 0):
            execute.WriteData(fid, prefix, obj = self, fieldname = 'subglacial_discharge', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        else:
            ## Write Integer fields
            fieldnames = ['sd_num_breaks', 'sd_num_params', 'sd_ar_order', 'sd_ma_order']
            for field in fieldnames:
                execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Integer')
            
            ## Write DoubleMat fields
            execute.WriteData(fid, prefix, obj = self, fieldname = 'sd_arma_timestep', format = 'Double', scale = md.constants.yts)
            execute.WriteData(fid, prefix, name = 'md.frontalforcings.sd_polynomialparams', data = sdpolyParams_scaled_2d, format = 'DoubleMat')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'sd_arlag_coefs',format = 'DoubleMat', yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname ='sd_malag_coefs', format = 'DoubleMat', yts = md.constants.yts)
            execute.WriteData(fid, prefix, name = 'md.frontalforcings.sd_datebreaks', data = sd_dbreaks, format = 'DoubleMat',scale = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'sd_monthlyfrac',format = 'DoubleMat', yts = md.constants.yts)