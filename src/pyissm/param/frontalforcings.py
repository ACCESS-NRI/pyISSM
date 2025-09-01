import numpy as np
from . import param_utils
from . import class_registry
from .. import execute

## ------------------------------------------------------
## frontalforcings.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    """
    Default frontalforcings parameters class for ISSM.

    This class encapsulates the default parameters for frontal forcings in the ISSM (Ice Sheet System Model) framework.
    It defines the main frontal forcing-related parameters.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    meltingrate : ndarray, default=np.nan
        Melting rate at given location [m/a].
    ablationrate : ndarray, default=np.nan
        Frontal ablation rate at given location [m/a], it contains both calving and melting.

    Methods
    -------
    __init__(self, other=None)
        Initializes the frontalforcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the frontalforcings parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.frontalforcings = pyissm.param.frontalforcings.default()
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

        s += '{}\n'.format(param_utils.fielddisplay(self, 'meltingrate', 'melting rate at given location [m/a]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ablationrate', 'frontal ablation rate at given location [m/a], it contains both calving and melting'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - frontalforcings.default Class'
        return s
    
    # Marshall method for saving the frontalforcings.default parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [frontalforcings.default] parameters to a binary file.

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
    Rignot frontalforcings parameters class for ISSM.

    This class encapsulates the parameters for frontal forcings based on the Rignot methodology in the ISSM (Ice Sheet System Model) framework.
    It defines the main frontal forcing-related parameters specific to the Rignot approach.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    basin_id : ndarray, default=np.nan
        Basin ID for elements.
    num_basins : int, default=0
        Number of basins.
    subglacial_discharge : ndarray, default=np.nan
        Sum of subglacial discharge for each basin [m/d].
    thermalforcing : ndarray, default=np.nan
        Thermal forcing [°C].

    Methods
    -------
    __init__(self, other=None)
        Initializes the Rignot frontalforcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the Rignot frontalforcings parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.frontalforcings = pyissm.param.frontalforcings.rignot()
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

        s += '{}\n'.format(param_utils.fielddisplay(self, 'basin_id', 'basin ID for elements'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'num_basins', 'number of basins'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'subglacial_discharge', 'sum of subglacial discharge for each basin [m/d]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'thermalforcing', 'thermal forcing [∘C]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - frontalforcings.rignot Class'
        return s
    
    # Marshall method for saving the frontalforcings.rignot parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [frontalforcings.rignot] parameters to a binary file.

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
        execute.WriteData(fid, prefix, name = 'md.frontalforcings.parameterization', data = 2, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, name = 'md.frontalforcings.basin_id', data = self.basin_id - 0, format = 'IntMat', mattype = 2) # 0-indexed
        execute.WriteData(fid, prefix, obj = self, fieldname = 'num_basins', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'subglacial_discharge', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'thermalforcing', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

## ------------------------------------------------------
## frontalforcings.rignotarma
## ------------------------------------------------------
@class_registry.register_class
class rignotarma(class_registry.manage_state):
    """
    RignotARMA frontalforcings parameters class for ISSM.

    This class encapsulates the parameters for frontal forcings based on the Rignot methodology with ARMA (AutoRegressive Moving Average) modeling in the ISSM (Ice Sheet System Model) framework.
    It defines the main frontal forcing-related parameters specific to the RignotARMA approach, including polynomial trends, breakpoints, ARMA coefficients, and subglacial discharge modeling.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    num_basins : int, default=0
        Number of different basins.
    num_params : int, default=0
        Number of different parameters in the piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.).
    num_breaks : int, default=0
        Number of different breakpoints in the piecewise-polynomial (separating num_breaks+1 periods).
    polynomialparams : ndarray, default=np.nan
        Coefficients for the polynomial (const, trend, quadratic, etc.), dim1 for basins, dim2 for periods, dim3 for orders.
    datebreaks : ndarray, default=np.nan
        Dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr].
    ar_order : int, default=0
        Order of the autoregressive model.
    ma_order : int, default=0
        Order of the moving-average model.
    arma_timestep : int, default=0
        Time resolution of the ARMA model [yr].
    arlag_coefs : ndarray, default=np.nan
        Basin-specific vectors of AR lag coefficients.
    malag_coefs : ndarray, default=np.nan
        Basin-specific vectors of MA lag coefficients.
    monthlyvals_intercepts : ndarray, default=np.nan
        Monthly intercept values for each basin.
    monthlyvals_trends : ndarray, default=np.nan
        Monthly trend values for each basin.
    monthlyvals_numbreaks : int, default=0
        Number of breakpoints for monthly values.
    monthlyvals_datebreaks : ndarray, default=np.nan
        Dates at which the monthly value breakpoints occur.
    basin_id : ndarray, default=np.nan
        Basin number assigned to each element.
    subglacial_discharge : ndarray, default=np.nan
        Sum of subglacial discharge for each basin [m/d].
    isdischargearma : int, default=0
        Whether an ARMA model is also used for the subglacial discharge (if 0: subglacial_discharge is used, if 1: sd_ parameters are used).
    sd_ar_order : int, default=0
        Order of the subglacial discharge autoregressive model.
    sd_ma_order : int, default=0
        Order of the subglacial discharge moving-average model.
    sd_arma_timestep : int, default=0
        Time resolution of the subglacial discharge ARMA model [yr].
    sd_arlag_coefs : ndarray, default=np.nan
        Basin-specific vectors of AR lag coefficients for subglacial discharge.
    sd_malag_coefs : ndarray, default=np.nan
        Basin-specific vectors of MA lag coefficients for subglacial discharge.
    sd_monthlyfrac : ndarray, default=np.nan
        Basin-specific vectors of 12 values with fraction of the annual discharge occurring every month.
    sd_num_breaks : int, default=0
        Number of different breakpoints in the subglacial discharge piecewise-polynomial (separating sd_num_breaks+1 periods).
    sd_num_params : int, default=0
        Number of different parameters in the subglacial discharge piecewise-polynomial.
    sd_polynomialparams : ndarray, default=np.nan
        Coefficients for the subglacial discharge polynomial (const, trend, quadratic, etc.).
    sd_datebreaks : ndarray, default=np.nan
        Dates at which the breakpoints in the subglacial discharge piecewise polynomial occur (1 row per basin) [yr].

    Methods
    -------
    __init__(self, other=None)
        Initializes the RignotARMA frontalforcings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the RignotARMA frontalforcings parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.frontalforcings = pyissm.param.frontalforcings.rignotarma()
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

        s += '{}\n'.format(param_utils.fielddisplay(self, 'num_basins', 'number of different basins [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'basin_id', 'basin number assigned to each element [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'num_breaks', 'number of different breakpoints in the piecewise-polynomial (separating num_breaks+1 periods)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'num_params', 'number of different parameters in the piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'polynomialparams', 'coefficients for the polynomial (const,trend,quadratic,etc.),dim1 for basins,dim2 for periods,dim3 for orders, ex: polyparams=cat(num_params,intercepts,trendlinearcoefs,trendquadraticcoefs)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'datebreaks', 'dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ar_order', 'order of the autoregressive model [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ma_order', 'order of the moving-average model [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'arma_timestep', 'time resolution of the ARMA model [yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'arlag_coefs', 'basin-specific vectors of AR lag coefficients [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'malag_coefs', 'basin-specific vectors of MA lag coefficients [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isdischargearma','whether an ARMA model is also used for the subglacial discharge (if 0: subglacial_discharge is used, if 1: sd_ parameters are used)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'subglacial_discharge', 'sum of subglacial discharge for each basin [m/d]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sd_ar_order','order of the subglacial discharge autoregressive model [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sd_ma_order','order of the subglacial discharge moving-average model [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sd_arma_timestep','time resolution of the subglacial discharge autoregressive model [yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sd_arlag_coefs','basin-specific vectors of AR lag coefficients for subglacial discharge [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sd_malag_coefs','basin-specific vectors of MA lag coefficients for subglacial discharge [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sd_monthlyfrac','basin-specific vectors of 12 values with fraction of the annual discharge occuring every month [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sd_num_params','number of different parameters in the subglacial discharge piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sd_num_breaks','number of different breakpoints in the subglacial discharge piecewise-polynomial (separating sd_num_breaks+1 periods)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sd_datebreaks','dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sd_polynomialparams','coefficients for the sd_polynomial (const,trend,quadratic,etc.),dim1 for basins,dim2 for periods,dim3 for orders'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - frontalforcings.rignotarma Class'
        return s
    
    # Marshall method for saving the frontalforcings.rignotarma parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [frontalforcings.rignotarma] parameters to a binary file.

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
        ## NOTE: Scaling logic here taken from $ISSM_DIR/src/m/classes/frontalforcingsrignotarma.py
        polyParams_scaled = np.copy(md.frontalforcings.polynomialparams)
        polyParams_scaled_2d = np.zeros((md.frontalforcings.num_basins, md.frontalforcings.num_breaks + 1 * md.frontalforcings.num_params))

        if(md.frontalforcings.num_params > 1):
            # Case 3D #
            if(md.frontalforcings.num_basins > 1 and md.frontalforcings.num_breaks + 1 > 1):
                for ii in range(md.frontalforcings.num_params):
                    polyParams_scaled[:,:,ii] = polyParams_scaled[:, :, ii] * (1 / md.constants.yts) ** ii
                # Fit in 2D array #
                for ii in range(md.frontalforcings.num_params):
                    polyParams_scaled_2d[:, ii * (md.frontalforcings.num_breaks + 1):(ii + 1) * (md.frontalforcings.num_breaks + 1)] = 1 * polyParams_scaled[:, :, ii]
            # Case 2D and higher-order params at increasing row index #
            elif(md.frontalforcings.num_basins == 1):
                for ii in range(md.frontalforcings.num_params):
                    polyParams_scaled[ii, :] = polyParams_scaled[ii, :] * (1 / md.constants.yts) ** ii
                # Fit in row array #
                for ii in range(md.frontalforcings.num_params):
                    polyParams_scaled_2d[0, ii * (md.frontalforcings.num_breaks + 1):(ii + 1) * (md.frontalforcings.num_breaks + 1)] = 1 * polyParams_scaled[ii, :]
            # Case 2D and higher-order params at increasing column index #
            elif(md.frontalforcings.num_breaks + 1 == 1):
                for ii in range(md.frontalforcings.num_params):
                    polyParams_scaled[:, ii] = polyParams_scaled[:, ii] * (1 / md.constants.yts) ** ii
                # 2D array is already in correct format #
                polyParams_scaled_2d = np.copy(polyParams_scaled)
        else:
            # 2D array is already in correct format and no need for scaling #
            polyParams_scaled_2d = np.copy(polyParams_scaled)
        if(md.frontalforcings.num_breaks + 1 == 1):
            dbreaks = np.zeros((md.frontalforcings.num_basins, 1))
        else:
            dbreaks = np.copy(md.frontalforcings.datebreaks)

        ### Deal with monthly effects ###
        if(np.any(np.isnan(md.frontalforcings.monthlyvals_intercepts))):
            interceptsM = np.zeros((md.frontalforcings.num_basins,12)) #monthly intercepts not provided, set to 0
            trendsM     = np.zeros((md.frontalforcings.num_basins,12)) #set monthly trends also to 0
        else:
            interceptsM3d = md.frontalforcings.monthlyvals_intercepts
            if(np.any(np.isnan(md.frontalforcings.monthlyvals_trends))):
                trendsM3d = 0*interceptsM3d #monthly trends not provided, set to 0
            else:
                trendsM3d = md.frontalforcings.monthlyvals_trends
        # Create 2D arrays from 3D arrays if needed #
        if(md.frontalforcings.monthlyvals_numbreaks + 1 > 1 and np.all(np.isnan(md.frontalforcings.monthlyvals_intercepts))==False):
            interceptsM = np.zeros((md.frontalforcings.num_basins, 12 * md.frontalforcings.monthlyvals_numbreaks + 1)) 
            trendsM     = np.zeros((md.frontalforcings.num_basins, 12 * md.frontalforcings.monthlyvals_numbreaks + 1))
            for ii in range(md.frontalforcings.monthlyvals_numbreaks + 1):
                interceptsM[:, ii * 12 : (ii + 1) * 12] = 1 * interceptsM3d[:,:,ii]
                trendsM[:, ii * 12 : (ii + 1) * 12] = 1 * trendsM3d[:,:,ii]
        elif(md.frontalforcings.monthlyvals_numbreaks + 1 == 1 and np.all(np.isnan(md.frontalforcings.monthlyvals_intercepts)) == False):
            interceptsM = 1 * interceptsM3d
            trendsM     = 1 * trendsM3d
        if(md.frontalforcings.monthlyvals_numbreaks + 1 == 1):
            dMbreaks = np.zeros((md.frontalforcings.num_basins, 1))
        else:
            dMbreaks = np.copy(md.frontalforcings.monthlyvals_datebreaks)

        ### Deal with the subglacial discharge polynomial ###
        if(self.isdischargearma):
            sdnprm  = md.frontalforcings.sd_num_params
            sdnper  = md.frontalforcings.sd_num_breaks+1
            sdpolyParams_scaled   = np.copy(md.frontalforcings.sd_polynomialparams)
            sdpolyParams_scaled_2d = np.zeros((md.frontalforcings.num_basins, md.frontalforcings.sd_num_breaks + 1 * md.frontalforcings.sd_num_params))
            if(md.frontalforcings.sd_num_params > 1):
                # Case 3D #
                if(nbas>1 and md.frontalforcings.sd_num_breaks + 1 > 1):
                    for ii in range(md.frontalforcings.sd_num_params):
                        sdpolyParams_scaled[:, :, ii] = sdpolyParams_scaled[:, :, ii] * (1 / md.constants.yts) ** ii
                    # Fit in 2D array #
                    for ii in range(md.frontalforcings.sd_num_params):
                        sdpolyParams_scaled_2d[:, ii * md.frontalforcings.sd_num_breaks + 1 : (ii + 1) * md.frontalforcings.sd_num_breaks + 1] = 1 * sdpolyParams_scaled[:, :, ii]
                # Case 2D and higher-order params at increasing row index #
                elif(nbas == 1):
                    for ii in range(md.frontalforcings.sd_num_params):
                        sdpolyParams_scaled[ii, :] = sdpolyParams_scaled[ii, :] * (1 / md.constants.yts) ** ii
                    # Fit in row array #
                    for ii in range(nprm):
                        sdpolyParams_scaled_2d[0, ii * md.frontalforcings.sd_num_breaks + 1 : (ii + 1) * md.frontalforcings.sd_num_breaks + 1] = 1 * sdpolyParams_scaled[ii, :]
                # Case 2D and higher-order params at incrasing column index #
                elif(md.frontalforcings.sd_num_breaks + 1 == 1):
                    for ii in range(md.frontalforcings.sd_num_params):
                        sdpolyParams_scaled[:, ii] = sdpolyParams_scaled[:, ii] * (1 / md.constants.yts) ** ii
                    # 2D array is already in correct format #
                    sdpolyParams_scaled_2d = np.copy(sdpolyParams_scaled)
            else:
                # 2D array is already in correct format and no need for scaling #
                sdpolyParams_scaled_2d = np.copy(sdpolyParams_scaled)
            if(md.frontalforcings.sd_num_breaks + 1 == 1):
                sd_dbreaks = np.zeros((md.frontalforcings.num_basins, 1))
            else:
                sd_dbreaks = np.copy(md.frontalforcings.sd_datebreaks)


        ## Write headers to file
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.frontalforcings.parameterization', data = 3, format = 'Integer')

        ## Write Integer fields
        fieldnames = ['num_basins', 'num_breaks', 'num_params', 'ar_order', 'ma_order', 'monthlyvals_numbreaks']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Integer')

        ## Write DoubleMat fields
        execute.WriteData(fid, prefix, name = 'md.frontalforcings.polynomialparams', data = polyparams2dScaled, format = 'DoubleMat')
        execute.WriteData(fid, prefix, 'object', self, fieldname = 'arlag_coefs', format = 'DoubleMat', yts = md.constants.yts)
        execute.WriteData(fid, prefix, 'object', self, fieldname = 'malag_coefs', format = 'DoubleMat', yts = md.constants.yts)
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
                WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname',field,'format','Integer')
            
            ## Write DoubleMat fields
            execute.WriteData(fid, prefix, obj = self, fieldname = 'sd_arma_timestep', format = 'Double', scale = md.constants.yts)
            execute.WriteData(fid, prefix, name = 'md.frontalforcings.sd_polynomialparams', data = sdpolyParams_scaled_2d, format = 'DoubleMat')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'sd_arlag_coefs',format = 'DoubleMat', yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname ='sd_malag_coefs', format = 'DoubleMat', yts = md.constants.yts)
            execute.WriteData(fid, prefix, name = 'md.frontalforcings.sd_datebreaks', data = sd_dbreaks, format = 'DoubleMat',scale = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'sd_monthlyfrac',format = 'DoubleMat', yts = md.constants.yts)