import numpy as np
from . import param_utils
from . import class_registry
from .. import execute
from .. import model

## ------------------------------------------------------
## smb.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    """
    Default surface mass balance (SMB) parameters class for ISSM.

    This class encapsulates the default parameters for surface mass balance in the ISSM (Ice Sheet System Model) framework.
    It defines the main SMB-related parameters.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    mass_balance : ndarray, default=np.nan
        Surface mass balance [m/yr ice eq].
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    requested_outputs : list, default=['default']
        Additional outputs requested
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic (default), 1: Geometric, 2: Harmonic.

    Methods
    -------
    __init__(self, other=None)
        Initializes the SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the SMB parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.smb = pyissm.param.smb.default()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.mass_balance = np.nan
        self.steps_per_step = 1
        self.requested_outputs = ['default']
        self.averaging = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'mass_balance', 'surface mass balance [m/yr ice eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.default Class'
        return s
    
    # Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self,
                        md = None,
                        return_default_outputs = False):
        """
        Process requested outputs, expanding 'default' to appropriate outputs.

        Parameters
        ----------
        md : ISSM model object, optional
            Model object containing mesh information.
        return_default_outputs : bool, default=False
            Whether to also return the list of default outputs.
            
        Returns
        -------
        outputs : list
            List of output strings with 'default' expanded to actual output names.
        default_outputs : list, optional
            Returned only if `return_default_outputs=True`.
        """

        outputs = []

        ## Set default_outputs
        default_outputs = ['SmbMassBalance']

        ## Loop through all requested outputs
        for item in self.requested_outputs:
            
            ## Process default outputs
            if item == 'default':
                    outputs.extend(default_outputs)

            ## Append other requested outputs (not defaults)
            else:
                outputs.append(item)

        if return_default_outputs:
            return outputs, default_outputs
        return outputs
        
    # Marshall method for saving the smb.default parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [smb.default] parameters to a binary file.

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
        execute.WriteData(fid, prefix, name = 'md.smb.model', data = 1, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'mass_balance', format = 'DoubleMat', scale = 1. / md.constants.yts, mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, name = 'md.smb.requested_outputs', data = self.process_outputs(md), format = 'StringArray')

## ------------------------------------------------------
## smb.arma
## ------------------------------------------------------
@class_registry.register_class
class arma(class_registry.manage_state):
    """
    ARMA (AutoRegressive Moving Average) surface mass balance model for ISSM.

    This class implements an ARMA-based surface mass balance model that combines 
    autoregressive and moving average components with piecewise polynomial trends
    and elevation-dependent lapse rates for basin-specific SMB modeling.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    num_basins : int, default=0
        Number of different basins [unitless].
    num_params : int, default=0
        Number of different parameters in the piecewise-polynomial (1:intercept only, 
        2:with linear trend, 3:with quadratic trend, etc.).
    num_breaks : int, default=0
        Number of different breakpoints in the piecewise-polynomial (separating 
        num_breaks+1 periods).
    polynomialparams : ndarray, default=np.nan
        Coefficients for the polynomial (const,trend,quadratic,etc.), dim1 for basins,
        dim2 for periods, dim3 for orders.
    ar_order : float, default=0.0
        Order of the autoregressive model [unitless].
    ma_order : float, default=0.0
        Order of the moving-average model [unitless].
    arlag_coefs : ndarray, default=np.nan
        Basin-specific vectors of AR lag coefficients [unitless].
    malag_coefs : ndarray, default=np.nan
        Basin-specific vectors of MA lag coefficients [unitless].
    datebreaks : ndarray, default=np.nan
        Dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr].
    basin_id : ndarray, default=np.nan
        Basin number assigned to each element [unitless].
    lapserates : ndarray, default=np.nan
        Basin-specific SMB lapse rates applied in each elevation bin, 1 row per basin,
        1 column per bin, dimension 3 can be of size 12 to prescribe monthly varying 
        values [m ice eq yr^-1 m^-1].
    elevationbins : ndarray, default=np.nan
        Basin-specific separations between elevation bins, 1 row per basin, 1 column 
        per limit between bins [m].
    refelevation : ndarray, default=np.nan
        Basin-specific reference elevations at which SMB is calculated [m].
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : list, default=['default']
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the ARMA SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the ARMA SMB parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.smb = pyissm.param.smb.arma()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.num_basins = 0
        self.num_params = 0
        self.num_breaks = 0
        self.polynomialparams = np.nan
        self.ar_order = 0.0
        self.ma_order = 0.0
        self.arlag_coefs = np.nan
        self.malag_coefs = np.nan
        self.polynomialparams = np.nan
        self.datebreaks = np.nan
        self.basin_id = np.nan
        self.lapserates = np.nan
        self.elevationbins = np.nan
        self.refelevation = np.nan
        self.datebreaks = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = ['default']

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

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
        s += '{}\n'.format(param_utils.fielddisplay(self, 'lapserates', 'basin-specific SMB lapse rates applied in each elevation bin, 1 row per basin, 1 column per bin, dimension 3 can be of size 12 to prescribe monthly varying values [m ice eq yr^-1 m^-1] (default: no lapse rate)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'elevationbins', 'basin-specific separations between elevation bins, 1 row per basin, 1 column per limit between bins, dimension 3 can be of size 12 to prescribe monthly varying values [m] (default: no basin separation)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'refelevation', 'basin-specific reference elevations at which SMB is calculated, and from which SMB is downscaled using lapserates (default: basin mean elevation) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.arma Class'
        return s

    # Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self,
                        md = None,
                        return_default_outputs = False):
        """
        Process requested outputs, expanding 'default' to appropriate outputs.

        Parameters
        ----------
        md : ISSM model object, optional
            Model object containing mesh information.
        return_default_outputs : bool, default=False
            Whether to also return the list of default outputs.
            
        Returns
        -------
        outputs : list
            List of output strings with 'default' expanded to actual output names.
        default_outputs : list, optional
            Returned only if `return_default_outputs=True`.
        """

        outputs = []

        ## Set default_outputs
        default_outputs = ['SmbMassBalance']

        ## Loop through all requested outputs
        for item in self.requested_outputs:
            
            ## Process default outputs
            if item == 'default':
                    outputs.extend(default_outputs)

            ## Append other requested outputs (not defaults)
            else:
                outputs.append(item)

        if return_default_outputs:
            return outputs, default_outputs
        return outputs
        
    # Marshall method for saving the smb.arma parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [smb.arma] parameters to a binary file.

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

        ## Scale parameters & set elevation bins
        ## NOTE: Taken from $ISSM_DIR/src/m/classes/SMBarma.py
        if(np.any(np.isnan(self.lapserates))):
            temp_lapse_rates = np.zeros((self.num_basins, 2, 12))
            print('      smb.lapserates not specified: set to 0')
            temp_elevation_bins = np.zeros((self.num_basins, 1, 12)) # Dummy elevation bins
            nbins    = 2
            ntmlapse = 12
        else:
            if len(np.shape(self.lapserates)) == 1:
                nbins    = 1
                ntmlapse = 1
            elif len(np.shape(self.lapserates)) == 2:
                nbins    = np.shape(self.lapserates)[1]
                ntmlapse = 1
            elif len(np.shape(self.lapserates)) == 3:
                nbins    = np.shape(self.lapserates)[1]
                ntmlapse = np.shape(self.lapserates)[2]
            temp_lapse_rates    = np.reshape(self.lapserates,[self.num_basins, nbins, ntmlapse])
            temp_elevation_bins = np.reshape(self.elevationbins, [self.num_basins, max(1, nbins - 1), ntmlapse])
        temp_ref_elevation  = np.copy(self.refelevation)
        
        # Scale the parameters
        polyParams_scaled   = np.copy(self.polynomialparams)
        polyParams_scaled_2d = np.zeros((self.num_basins, self.num_breaks + 1 * self.num_params))
        if self.num_params > 1:
            # Case 3D
            if self.num_basins > 1 and self.num_breaks + 1 > 1:
                for ii in range(self.num_params):
                    polyParams_scaled[:, :, ii] = polyParams_scaled[:, :, ii] * (1 / md.constants.yts) ** (ii + 1)
                # Fit in 2D array
                for ii in range(self.num_params):
                    polyParams_scaled_2d[:, ii * self.num_breaks + 1:(ii + 1) * self.num_breaks + 1] = 1 * polyParams_scaled[:, :, ii]
            # Case 2D and higher-order params at increasing row index
            elif self.num_basins == 1:
                for ii in range(self.num_params):
                    polyParams_scaled[ii, :] = polyParams_scaled[ii, :] * (1 / md.constants.yts) ** (ii + 1)
                # Fit in row array
                for ii in range(self.num_params):
                    polyParams_scaled_2d[0, ii * self.num_breaks + 1:(ii + 1) * self.num_breaks + 1] = 1 * polyParams_scaled[ii, :]
            # Case 2D and higher-order params at increasing column index
            elif self.num_breaks + 1 == 1:
                for ii in range(self.num_params):
                    polyParams_scaled[:, ii] = polyParams_scaled[:, ii] * (1 / md.constants.yts) ** (ii + 1)
                # 2D array is already in correct format
                polyParams_scaled_2d = np.copy(polyParams_scaled)
        else:
            polyParams_scaled   = polyParams_scaled * (1 / md.constants.yts)
            # 2D array is already in correct format
            polyParams_scaled_2d = np.copy(polyParams_scaled)

        if self.num_breaks + 1 == 1:
            dbreaks = np.zeros((self.num_basins, 1))
        else:
            dbreaks = np.copy(self.datebreaks)

        if ntmlapse == 1:
            temp_lapse_rates    = np.repeat(temp_lapse_rates, 12, axis = 2)
            temp_elevation_bins = np.repeat(temp_elevation_bins, 12, axis = 2)
        if np.any(np.isnan(self.refelevation)):
            temp_ref_elevation = np.zeros((self.num_basins)).reshape(1, self.num_basins)
            areas = model.mesh.get_element_areas_volumes(md.mesh.elements, md.mesh.x, md.mesh.y)
            for ii, bid in enumerate(np.unique(self.basin_id)):
                indices = np.where(self.basin_id == bid)[0]
                elemsh = np.zeros((len(indices)))
                for jj in range(len(indices)):
                    elemsh[jj] = np.mean(md.geometry.surface[md.mesh.elements[indices[jj], :] - 1])
                temp_ref_elevation[0, ii] = np.sum(areas[indices] * elemsh) / np.sum(areas[indices])
            if(np.any(temp_lapse_rates != 0)):
                print('      smb.refelevation not specified: Reference elevations set to mean surface elevation of basins')
        nbins = np.shape(temp_lapse_rates)[1]
        temp_lapse_rates_2d    = np.zeros((self.num_basins, nbins * 12))
        temp_elevation_bins_2d = np.zeros((self.num_basins, max(12, (nbins - 1) * 12)))
        for ii in range(12):
            temp_lapse_rates_2d[:, ii * nbins:(ii + 1) * nbins] = temp_lapse_rates[:, :, ii]
            temp_elevation_bins_2d[:, ii * (nbins - 1):(ii + 1) * (nbins - 1)] = temp_elevation_bins[:, :, ii]

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.smb.model', data = 13, format = 'Integer')

        ## Write Integer fields
        fieldnames = ['num_basins', 'num_breaks', 'num_params', 'ar_order', 'ma_order', 'steps_per_step', 'averaging']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.smb.num_bins', data = nbins, format = 'Integer')

        ## Write DoubleMat fields
        execute.WriteData(fid, prefix, name = 'md.smb.polynomialparams', data = polyParams_scaled_2d,  format = 'DoubleMat')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'arlag_coefs', format = 'DoubleMat', yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'malag_coefs', format = 'DoubleMat', yts = md.constants.yts)
        execute.WriteData(fid, prefix, name = 'md.smb.datebreaks', data = dbreaks, format = 'DoubleMat', yts = md.constants.yts)
        execute.WriteData(fid, prefix, name = 'smb.smb.lapserates', data = temp_lapse_rates_2d, format = 'DoubleMat', scale = 1. / md.constants.yts, yts = md.constants.yts)
        execute.WriteData(fid, prefix, name = 'md.smb.elevationbins', data = temp_elevation_bins_2d, format = 'DoubleMat')
        execute.WriteData(fid, prefix, name = 'md.smb.refelevation', data = temp_ref_elevation, format = 'DoubleMat')

        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'arma_timestep', format = 'Double', scale = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, name = 'md.smb.basin_id', data = self.basin_id - 1, format = 'IntMat', mattype = 2)  # 0-indexed
        execute.WriteData(fid, prefix, name = 'md.smb.requested_outputs', data = self.process_outputs(md), format = 'StringArray')

## ------------------------------------------------------
## smb.components
## ------------------------------------------------------
@class_registry.register_class
class components(class_registry.manage_state):
    """
    Component-based surface mass balance model for ISSM.

    This class implements a component-based SMB model where the surface mass balance
    is calculated as SMB = accumulation - runoff - evaporation. Each component can
    be specified independently.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    accumulation : ndarray, default=np.nan
        Accumulated snow [m/yr ice eq].
    runoff : ndarray, default=np.nan
        Amount of ice melt lost from the ice column [m/yr ice eq].
    evaporation : ndarray, default=np.nan
        Amount of ice lost to evaporative processes [m/yr ice eq].
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : list, default=['default']
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the component SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the component SMB parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Notes
    -----
    The surface mass balance is computed as:
    SMB = accumulation - runoff - evaporation

    Examples
    --------
    md.smb = pyissm.param.smb.components()
    md.smb.accumulation = accumulation_data
    md.smb.runoff = runoff_data
    md.smb.evaporation = evaporation_data
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.accumulation = np.nan
        self.runoff = np.nan
        self.evaporation = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = ['default']

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters (SMB=accumulation-runoff-evaporation) :\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'accumulation', 'accumulated snow [m/yr ice eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'runoff', 'amount of ice melt lost from the ice column [m/yr ice eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'evaporation', 'mount of ice lost to evaporative processes [m/yr ice eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.components Class'
        return s
    
    # Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self,
                        md = None,
                        return_default_outputs = False):
        """
        Process requested outputs, expanding 'default' to appropriate outputs.

        Parameters
        ----------
        md : ISSM model object, optional
            Model object containing mesh information.
        return_default_outputs : bool, default=False
            Whether to also return the list of default outputs.
            
        Returns
        -------
        outputs : list
            List of output strings with 'default' expanded to actual output names.
        default_outputs : list, optional
            Returned only if `return_default_outputs=True`.
        """

        outputs = []

        ## Set default_outputs
        default_outputs = ['SmbMassBalance']

        ## Loop through all requested outputs
        for item in self.requested_outputs:
            
            ## Process default outputs
            if item == 'default':
                    outputs.extend(default_outputs)

            ## Append other requested outputs (not defaults)
            else:
                outputs.append(item)

        if return_default_outputs:
            return outputs, default_outputs
        return outputs
        
    # Marshall method for saving the smb.components parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [smb.components] parameters to a binary file.

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
        execute.WriteData(fid, prefix, name = 'md.smb.model', data = 2, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'accumulation', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'runoff', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'evaporation', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.smb.requested_outputs', data = self.process_outputs(md), format = 'StringArray')


## ------------------------------------------------------
## smb.d18opdd
## ------------------------------------------------------
@class_registry.register_class
class d18opdd(class_registry.manage_state):
    """
    Delta-18-O driven positive degree day surface mass balance model for ISSM.

    This class implements a positive degree day (PDD) SMB model driven by delta-18-O 
    isotope data for paleoclimate applications. It includes temperature and precipitation
    scaling based on isotope ratios and elevation-dependent corrections.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    desfac : float, default=0.5
        Desertification elevation factor (between 0 and 1) [m].
    s0p : ndarray, default=np.nan
        Elevation from precipitation source (between 0 and a few 1000s m) [m].
    s0t : ndarray, default=np.nan
        Elevation from temperature source (between 0 and a few 1000s m) [m].
    rlaps : float, default=6.5
        Present day lapse rate [degree/km].
    rlapslgm : float, default=6.5
        LGM lapse rate [degree/km].
    dpermil : float, default=2.4
        Degree per mil, required if d18opd is activated.
    f : float, default=0.169
        Precipitation/temperature scaling factor, required if d18opd is activated.
    Tdiff : ndarray, default=np.nan
        Temperature difference field.
    sealev : ndarray, default=np.nan
        Sea level data.
    ismungsm : int, default=0
        Is mungsm parametrisation activated (0 or 1).
    isd18opd : int, default=1
        Is delta18o parametrisation from present day temperature and precipitation activated (0 or 1).
    issetpddfac : int, default=0
        Is user passing in defined PDD factors (0 or 1).
    istemperaturescaled : int, default=1
        Is temperature scaled to delta18o value (0 or 1).
    isprecipscaled : int, default=1
        Is precipitation scaled to delta18o value (0 or 1).
    delta18o : ndarray, default=np.nan
        Delta-18-O values [per mil].
    delta18o_surface : ndarray, default=np.nan
        Surface delta-18-O values.
    temperatures_presentday : ndarray, default=np.nan
        Monthly present day surface temperatures [K].
    precipitations_presentday : ndarray, default=np.nan
        Monthly surface precipitation [m/yr water eq].
    temperatures_reconstructed : ndarray, default=np.nan
        Monthly historical surface temperatures [K].
    precipitations_reconstructed : ndarray, default=np.nan
        Monthly historical precipitation [m/yr water eq].
    pddfac_snow : ndarray, default=np.nan
        PDD factor for snow [mm ice equiv/day/degree C].
    pddfac_ice : ndarray, default=np.nan
        PDD factor for ice [mm ice equiv/day/degree C].
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : list, default=['default']
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the delta-18-O PDD SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the delta-18-O PDD SMB parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.smb = pyissm.param.smb.d18opdd()
    md.smb.delta18o = delta18o_data
    md.smb.temperatures_presentday = temp_data
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.desfac = 0.5
        self.s0p = np.nan
        self.s0t = np.nan
        self.rlaps = 6.5
        self.rlapslgm = 6.5
        self.dpermil = 2.4
        self.f = 0.169
        self.Tdiff = np.nan
        self.sealev = np.nan
        self.ismungsm = 0
        self.isd18opd = 1
        self.issetpddfac = 0
        self.istemperaturescaled = 1
        self.isprecipscaled = 1
        self.delta18o = np.nan
        self.delta18o_surface = np.nan
        self.temperatures_presentday = np.nan
        self.precipitations_presentday = np.nan
        self.temperatures_reconstructed = np.nan
        self.precipitations_reconstructed = np.nan
        self.pddfac_snow = np.nan
        self.pddfac_ice = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'isd18opd', 'is delta18o parametrisation from present day temperature and precipitation activated (0 or 1, default is 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'issetpddfac', 'is user passing in defined pdd factors (0 or 1, default is 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'desfac', 'desertification elevation factor (between 0 and 1, default is 0.5) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 's0p', 'should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 's0t', 'should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rlaps', 'present day lapse rate [degree/km]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'temperatures_presentday', 'monthly present day surface temperatures [K], required if delta18o/mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'precipitations_presentday', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'istemperaturescaled', 'if delta18o parametrisation from present day temperature and precipitation is activated, is temperature scaled to delta18o value? (0 or 1, default is 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isprecipscaled', 'if delta18o parametrisation from present day temperature and precipitation is activated, is precipitation scaled to delta18o value? (0 or 1, default is 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'temperatures_reconstructed', 'monthly historical surface temperatures [K], required if delta18o/mungsm/d18opd is activated and istemperaturescaled is not activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'precipitations_reconstructed', 'monthly historical precipitation [m/yr water eq], required if delta18o/mungsm/d18opd is activated and isprecipscaled is not activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'delta18o', 'delta18o [per mil], required if pdd is activated and delta18o activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'dpermil', 'degree per mil, required if d18opd is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'f', 'precip/temperature scaling factor, required if d18opd is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'pddfac_snow', 'Pdd factor for snow for all the domain [mm ice equiv/day/degree C]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'pddfac_ice', 'Pdd factor for ice for all the domain [mm ice equiv/day/degree C]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.d18opdd Class'
        return s
    
    # Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self,
                        md = None,
                        return_default_outputs = False):
        """
        Process requested outputs, expanding 'default' to appropriate outputs.

        Parameters
        ----------
        md : ISSM model object, optional
            Model object containing mesh information.
        return_default_outputs : bool, default=False
            Whether to also return the list of default outputs.
            
        Returns
        -------
        outputs : list
            List of output strings with 'default' expanded to actual output names.
        default_outputs : list, optional
            Returned only if `return_default_outputs=True`.
        """

        outputs = []

        ## Set default_outputs
        default_outputs = ['SmbMassBalance']

        ## Loop through all requested outputs
        for item in self.requested_outputs:
            
            ## Process default outputs
            if item == 'default':
                    outputs.extend(default_outputs)

            ## Append other requested outputs (not defaults)
            else:
                outputs.append(item)

        if return_default_outputs:
            return outputs, default_outputs
        return outputs
        
    # Marshall method for saving the smb.d18opdd parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [smb.d18opdd] parameters to a binary file.

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
        execute.WriteData(fid, prefix, name = 'md.smb.model', data = 5, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'ismungsm', format = 'Boolean')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'isd18opd', format = 'Boolean')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'issetpddfac', format = 'Boolean')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'desfac', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 's0p', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 's0t', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rlaps', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rlapslgm', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'Tdiff', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'sealev', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.smb.requested_outputs', data = self.process_outputs(md), format = 'StringArray')

        ## Write conditional fields
        if self.isd18opd:
            execute.WriteData(fid, prefix, obj = self, fieldname = 'temperatures_presentday', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'precipitations_presentday', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'istemperaturescaled', format = 'Boolean')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'isprecipscaled', format = 'Boolean')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'delta18o', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'dpermil', format = 'Double')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'f', format = 'Double')

            if self.istemperaturescaled == 0:
                execute.WriteData(fid, prefix, obj = self, fieldname = 'temperatures_reconstructed', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

            if self.isprecipscaled == 0:
                execute.WriteData(fid, prefix, obj = self, fieldname = 'precipitations_reconstructed', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

        if self.issetpddfac:
            execute.WriteData(fid, prefix, obj = self, fieldname = 'pddfac_snow', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'pddfac_ice', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

## ------------------------------------------------------
## smb.gemb
## ------------------------------------------------------
@class_registry.register_class
class gemb(class_registry.manage_state):
    """
    GEMB (Greenland Energy and Mass Balance) surface mass balance model for ISSM.

    This class implements the GEMB model, a sophisticated physics-based surface mass 
    balance model that simulates snow/firn densification, grain evolution, albedo, 
    shortwave radiation penetration, thermal processes, melt, and accumulation. 
    Originally developed for Greenland but applicable to other ice sheets.

    Parameters
    ----------
    md : ISSM model object
        Model object containing mesh information. Required for initialization.
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    isgraingrowth : int, default=1
        Run grain growth module (1=on, 0=off).
    isalbedo : int, default=1
        Run albedo module (1=on, 0=off).
    isshortwave : int, default=1
        Run short wave radiation module (1=on, 0=off).
    isthermal : int, default=1
        Run thermal module (1=on, 0=off).
    isaccumulation : int, default=1
        Run accumulation module (1=on, 0=off).
    ismelt : int, default=1
        Run melting module (1=on, 0=off).
    isdensification : int, default=1
        Run densification module (1=on, 0=off).
    isturbulentflux : int, default=1
        Run turbulent heat fluxes module (1=on, 0=off).
    isconstrainsurfaceT : int, default=1
        Constrain surface temperatures to air temperature (1=on, 0=off).
    isdeltaLWup : int, default=0
        Apply bias to long wave upward radiation spatially (1=on, 0=off).
    ismappedforcing : int, default=0
        Use mapped forcing when grid doesn't match model mesh (1=on, 0=off).
    iscompressedforcing : int, default=0
        Compress input matrices when writing to binary (1=on, 0=off).
    Ta : ndarray, default=np.nan
        2-meter air temperature [K].
    V : ndarray, default=np.nan
        Wind speed [m/s].
    dswrf : ndarray, default=np.nan
        Downward shortwave radiation flux [W/m²].
    dlwrf : ndarray, default=np.nan
        Downward longwave radiation flux [W/m²].
    P : ndarray, default=np.nan
        Precipitation [mm w.e./m²].
    eAir : ndarray, default=np.nan
        Screen level vapor pressure [Pa].
    pAir : ndarray, default=np.nan
        Surface pressure [Pa].
    Tmean : ndarray, default=np.nan
        Mean annual temperature [K].
    Vmean : ndarray, default=10.0
        Mean annual wind speed [m/s].
    C : ndarray, default=np.nan
        Mean annual snow accumulation [kg/m²/yr].
    Tz : ndarray, default=np.nan
        Height above ground at which temperature was sampled [m].
    Vz : ndarray, default=np.nan
        Height above ground at which wind was sampled [m].
    zTop : ndarray, default=10.0
        Depth over which grid length is constant at snowpack top [m].
    dzTop : ndarray, default=0.05
        Initial top vertical grid spacing [m].
    dzMin : ndarray, default=dzTop/2
        Initial minimum allowable vertical grid spacing [m].
    zY : ndarray, default=1.025
        Grid stretching factor below top zone.
    zMax : ndarray, default=250.0
        Initial maximum model depth [m].
    zMin : ndarray, default=130.0
        Initial minimum model depth [m].
    aIdx : int, default=1
        Albedo calculation method (0-4).
    eIdx : int, default=1
        Emissivity calculation method (0-2).
    tcIdx : int, default=1
        Thermal conductivity method (1-2).
    swIdx : int, default=0
        Shortwave penetration method (0-1).
    denIdx : int, default=2
        Densification model (1-7).
    dsnowIdx : int, default=1
        Fresh snow density model (0-4).
    outputFreq : int, default=30
        Output frequency [days].
    InitDensityScaling : float, default=1.0
        Initial density scaling factor.
    ThermoDeltaTScaling : float, default=1/11.0
        Thermal diffusion timestep scaling factor.
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : list, default=[]
        Additional outputs requested.

    Methods
    -------
    __init__(self, md=None, other=None)
        Initializes the GEMB SMB parameters, requires model object for mesh information.
    __repr__(self)
        Returns a detailed string representation of the GEMB SMB parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Notes
    -----
    GEMB is a comprehensive physically-based surface mass balance model that includes:
    - Multi-layer snow/firn column with evolving properties
    - Grain size evolution and metamorphism
    - Surface albedo calculation with multiple parameterizations
    - Shortwave radiation penetration and absorption
    - Thermal diffusion and temperature evolution
    - Melt and refreeze processes
    - Snow densification using various empirical/physical models
    - Turbulent heat flux calculations

    The model requires detailed meteorological forcing (temperature, wind, radiation,
    precipitation, humidity, pressure) and simulates the evolution of snow/firn 
    properties over time.

    Examples
    --------
    md.smb = pyissm.param.smb.gemb(md)
    md.smb.Ta = temperature_forcing_data
    md.smb.P = precipitation_data
    md.smb.dswrf = shortwave_radiation_data
    """

    # Initialise with default parameters
    ## NOTE: md must be specified for mesh parameters.
    def __init__(self,
                 md = None,
                 other = None):
        self.isgraingrowth = 1
        self.isalbedo = 1
        self.isshortwave = 1
        self.isthermal = 1
        self.isaccumulation = 1
        self.ismelt = 1
        self.isdensification = 1
        self.isturbulentflux = 1
        self.isconstrainsurfaceT = 1
        self.isdeltaLWup = 0
        self.ismappedforcing = 0
        self.iscompressedforcing = 0
        self.Ta = np.nan
        self.V = np.nan
        self.dswrf = np.nan
        self.dlwrf = np.nan
        self.P = np.nan
        self.eAir = np.nan
        self.pAir = np.nan
        self.Tmean = np.nan
        self.Vmean = 10 * np.ones((md.mesh.numberofelements,))
        self.C = np.nan
        self.Tz = np.nan
        self.Vz = np.nan
        self.aValue = self.aSnow * np.ones(md.mesh.numberofelements,)
        self.teValue =  np.ones((md.mesh.numberofelements,))
        self.dulwrfValue = np.ones((md.mesh.numberofelements,))
        self.mappedforcingpoint = np.nan
        self.mappedforcingelevation = np.nan
        self.lapseTaValue = -0.006
        self.lapsedlwrfValue = -0.032
        self.Dzini = 0.05 * np.ones((md.mesh.numberofelements, 2))
        self.Dini = 910.0 * np.ones((md.mesh.numberofelements, 2))
        self.Reini = 2.5 * np.ones((md.mesh.numberofelements, 2))
        self.Gdnini = 0.0 * np.ones((md.mesh.numberofelements, 2))
        self.Gspini = 0.0 * np.ones((md.mesh.numberofelements, 2))
        self.ECini = 0.0 * np.ones((md.mesh.numberofelements, 2))
        self.Wini = 0.0 * np.ones((md.mesh.numberofelements, 2))
        self.Aini = 0.0 * np.ones((md.mesh.numberofelements, 2))
        self.Adiffini = 0.0 * np.ones((md.mesh.numberofelements, 2))
        self.Tini = 273.15 * np.ones((md.mesh.numberofelements, 2))
        self.Sizeini = 2 * np.ones((md.mesh.numberofelements, ))
        self.aIdx = 1
        self.eIdx = 1
        self.tcIdx = 1
        self.swIdx = 0
        self.denIdx = 2
        self.dsnowIdx = 1
        self.zTop = 10 * np.ones((md.mesh.numberofelements,))
        self.dzTop = 0.05 * np.ones((md.mesh.numberofelements,))
        self.dzMin = self.dzTop / 2
        self.zY = 1.025 * np.ones((md.mesh.numberofelements,))
        self.zMax = 250 * np.ones((md.mesh.numberofelements,))
        self.zMin = 130  * np.ones((md.mesh.numberofelements,))
        self.outputFreq = 30
        self.dswdiffrf = 0.0 * np.ones(md.mesh.numberofelements,)
        self.szaValue = 0.0 * np.ones(md.mesh.numberofelements,)
        self.cotValue = 0.0 * np.ones(md.mesh.numberofelements,)
        self.ccsnowValue = 0.0 * np.ones(md.mesh.numberofelements,)
        self.cciceValue = 0.0 * np.ones(md.mesh.numberofelements,)
        self.aSnow = 0.85
        self.aIce = 0.48
        self.cldFrac = 0.1
        self.t0wet = 15
        self.t0dry = 30
        self.K = 7
        self.adThresh = 1023
        self.teThresh = 10
        self.InitDensityScaling = 1.0
        self.ThermoDeltaTScaling = 1 / 11.0
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = []

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings for SMB GEMB model :\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isgraingrowth', 'run grain growth module (default true)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isalbedo', 'run albedo module (default true)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isshortwave', 'run short wave module (default true)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isthermal', 'run thermal module (default true)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isaccumulation', 'run accumulation module (default true)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ismelt', 'run melting  module (default true)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isdensification', 'run densification module (default true)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isturbulentflux', 'run turbulant heat fluxes module (default true)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isconstrainsurfaceT', 'constrain surface temperatures to air temperature, turn off EC and surface flux contribution to surface temperature change (default false)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isdeltaLWup', 'set to true to invoke a bias in the long wave upward spatially, specified by dulwrfValue (default false)'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'ismappedforcing','set to true if forcing grid does not match model mesh, mapping specified by mappedforcingpoint (default false)'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'iscompressedforcing','set to true to compress the input matrices when writing to binary (default false)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Ta', '2 m air temperature, in Kelvin'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'V', 'wind speed (m s-1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'dswrf', 'downward shortwave radiation flux [W/m^2]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'dswdiffrf', 'downward diffusive portion of shortwave radiation flux (default to 0) [W/m^2]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'dlwrf', 'downward longwave radiation flux [W/m^2]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'P', 'precipitation [mm w.e. / m^2]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'eAir', 'screen level vapor pressure [Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'pAir', 'surface pressure [Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Tmean', 'mean annual temperature [K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'C', 'mean annual snow accumulation [kg m-2 yr-1]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Vmean', 'mean annual temperature [m s-1] (default 10 m/s)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Tz', 'height above ground at which temperature (T) was sampled [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Vz', 'height above ground at which wind (V) eas sampled [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'zTop', 'depth over which grid length is constant at the top of the snopack (default 10) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'dzTop', 'initial top vertical grid spacing (default .05) [m] '))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'dzMin', 'initial min vertical allowable grid spacing (default dzMin/2) [m] '))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'zMax', 'initial max model depth (default is min(thickness, 500)) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'zMin', 'initial min model depth (default is min(thickness, 30)) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'zY', 'stretch grid cells bellow top_z by a [top_dz * y ^ (cells bellow top_z)]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'InitDensityScaling', ['initial scaling factor multiplying the density of ice', 'which describes the density of the snowpack.']))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ThermoDeltaTScaling', 'scaling factor to multiply the thermal diffusion timestep (delta t)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'outputFreq', 'output frequency in days (default is monthly, 30)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'adThresh', 'Apply aIdx method to all areas with densities below this value, or else apply direct input value from aValue, allowing albedo to be altered.'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'aIdx', ['method for calculating albedo and subsurface absorption (default is 1)',
            '0: direct input from aValue parameter',
            '1: effective grain radius [Gardner & Sharp, 2009]',
            '2: effective grain radius [Brun et al., 1992; LeFebre et al., 2003], with swIdx=1, SW penetration follows grain size in 3 spectral bands (Brun et al., 1992)',
            '3: density and cloud amount [Greuell & Konzelmann, 1994]',
            '4: exponential time decay & wetness [Bougamont & Bamber, 2005]']))

        s += '{}\n'.format(param_utils.fielddisplay(self, 'dulwrfValue', 'Specified bias to be applied to the outward long wave radiation at every element (W/m-2, +upward)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'teValue', 'Outward longwave radiation thermal emissivity forcing at every element (default in code is 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'teThresh', ['Apply eIdx method to all areas with effective grain radius above this value (mm),', 'or else apply direct input value from teValue, allowing emissivity to be altered.']))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'eIdx', ['method for calculating emissivity (default is 1)',
            '0: direct input from teValue parameter, no use of teThresh',
            '1: default value of 1, in areas with grain radius below teThresh',
            '2: default value of 1, in areas with grain radius below teThresh and areas of dry snow (not bare ice or wet) at the surface']))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'tcIdx', ['method for calculating thermal conductivity (default is 1)',
            '1: after Sturm et al, 1997',
            '2: after Calonne et al., 2011']))

        s += '{}\n'.format(param_utils.fielddisplay(self,'mappedforcingpoint','Mapping of which forcing point will map to each mesh element for ismappedforcing option (integer). Size number of elements.'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'mappedforcingelevation','The elevation of each mapped forcing location (m above sea level) for ismappedforcing option. Size number of forcing points.'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'lapseTaValue','Temperature lapse rate if forcing has different grid and should be remapped for ismappedforcing option. (Default value is -0.006 K m-1.)'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'lapsedlwrfValue','Longwave down lapse rate if forcing has different grid and should be remapped for ismappedforcing option. (Default value is -0.032 W m-2 m-1.)'))

        # Snow properties init
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Dzini', 'Initial cell depth when restart [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Dini', 'Initial snow density when restart [kg m-3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Reini', 'Initial grain size when restart [mm]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Gdnini', 'Initial grain dricity when restart [-]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Gspini', 'Initial grain sphericity when restart [-]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ECini', 'Initial evaporation/condensation when restart [kg m-2]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Wini', 'Initial snow water content when restart [kg m-2]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Aini', 'Initial albedo when restart [-]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Adiffini', 'Initial diffusive radiation albedo when restart (default to 1) [-]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Tini', 'Initial snow temperature when restart [K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Sizeini', 'Initial number of layers when restart [-]'))

        # Additional albedo parameters
        s += '{}\n'.format(param_utils.fielddisplay(self, 'aValue', 'Albedo forcing at every element'))
        if isinstance(self.aIdx, (list, type(np.array([1, 2])))) and (self.aIdx == [1, 2] or (1 in self.aIdx and 2 in self.aIdx)):
            s += '{}\n'.format(param_utils.fielddisplay(self, 'aSnow', 'new snow albedo (0.64 - 0.89)'))
            s += '{}\n'.format(param_utils.fielddisplay(self, 'aIce', 'albedo of ice (0.27-0.58)'))
            if self.aIdx == 1:
                s += '{}\n'.format(param_utils.fielddisplay(self,'szaValue','Solar Zenith Angle [degree]'))
                s += '{}\n'.format(param_utils.fielddisplay(self,'cotValue','Cloud Optical Thickness'))
                s += '{}\n'.format(param_utils.fielddisplay(self,'ccsnowValue','concentration of light absorbing carbon for snow [ppm1]'))
                s += '{}\n'.format(param_utils.fielddisplay(self,'cciceValue','concentration of light absorbing carbon for ice [ppm1]'))
        elif self.aIdx == 3:
            s += '{}\n'.format(param_utils.fielddisplay(self, 'cldFrac', 'average cloud amount'))
        elif self.aIdx == 4:
            s += '{}\n'.format(param_utils.fielddisplay(self, 't0wet', 'time scale for wet snow (15-21.9) [d]'))
            s += '{}\n'.format(param_utils.fielddisplay(self, 't0dry', 'warm snow timescale (30) [d]'))
            s += '{}\n'.format(param_utils.fielddisplay(self, 'K', 'time scale temperature coef. (7) [d]'))

        s += '{}\n'.format(param_utils.fielddisplay(self, 'swIdx', 'apply all SW to top grid cell (0) or allow SW to penetrate surface (1) [default 0, if swIdx=1 and aIdx=2 function of effective radius (Brun et al., 1992) or else dependent on snow density (taken from Bassford, 2002)]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'denIdx', ['densification model to use (default is 2):',
            '1 = emperical model of Herron and Langway (1980)',
            '2 = semi-emperical model of Anthern et al. (2010)',
            '3 = DO NOT USE: physical model from Appix B of Anthern et al. (2010)',
            '4 = DO NOT USE: emperical model of Li and Zwally (2004)',
            '5 = DO NOT USE: modified emperical model (4) by Helsen et al. (2008)',
            '6 = Antarctica semi-emperical model of Ligtenberg et al. (2011)',
            '7 = Greenland semi-emperical model of Kuipers Munneke et al. (2015)']))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'dsnowIdx', ['model for fresh snow accumulation density (default is 1):',
            '0 = Original GEMB value, 150 kg/m^3',
            '1 = Antarctica value of fresh snow density, 350 kg/m^3',
            '2 = Greenland value of fresh snow density, 315 kg/m^3, Fausto et al. (2018)',
            '3 = Antarctica model of Kaspers et al. (2004), Make sure to set Vmean accurately',
            '4 = Greenland model of Kuipers Munneke et al. (2015)']))

        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s
    
    # Define class string
    def __str__(self):
        s = 'ISSM - smb.gemb Class'
        return s

    # Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self,
                        md = None,
                        return_default_outputs = False):
        """
        Process requested outputs, expanding 'default' to appropriate outputs.

        Parameters
        ----------
        md : ISSM model object, optional
            Model object containing mesh information.
        return_default_outputs : bool, default=False
            Whether to also return the list of default outputs.
            
        Returns
        -------
        outputs : list
            List of output strings with 'default' expanded to actual output names.
        default_outputs : list, optional
            Returned only if `return_default_outputs=True`.
        """

        outputs = []

        ## Set default_outputs
        default_outputs = ['SmbMassBalance','SmbAccumulatedMassBalance']

        ## Loop through all requested outputs
        for item in self.requested_outputs:
            
            ## Process default outputs
            if item == 'default':
                    outputs.extend(default_outputs)

            ## Append other requested outputs (not defaults)
            else:
                outputs.append(item)

        if return_default_outputs:
            return outputs, default_outputs
        return outputs

    # Marshall method for saving the smb.gemb parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [smb.gemb] parameters to a binary file.

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
        execute.WriteData(fid, prefix, name = 'md.smb.model', data = 8, format = 'Integer')

        ## Write Boolean fields
        fieldnames = ['isgraingrowth', 'isalbedo', 'isshortwave', 'isthermal', 'isaccumulation',
                      'ismelt', 'isdensification', 'isturbulentflux', 'isconstrainsurfaceT',
                      'isdeltaLWup', 'ismappedforcing']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Boolean')

        ## Write conditional compressed forcing fields
        if self.iscompressedforcing:
            writetype='CompressedMat'
        else:
            writetype='DoubleMat'

        fieldnames = ['Ta', 'V', 'dswrf', 'dswdiffrf', 'dlwrf', 'P', 'eAir', 'pAir']
        for field in fieldnames:
            execute.WriteData(fid,prefix, obj = self, fieldname = field, format = writetype, mattype = 2, timeserieslength = np.shape(self.Ta)[0], yts = md.constants.yts)

        ## Write DoubleMat fields
        fieldnames = ['Tmean', 'C', 'Vmean', 'Tz', 'Vz', 'zTop', 'dzTop', 'dzMin', 'zY', 'zMax', 'zMin']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'DoubleMat', mattype = 2)

        ## Write Integer fields
        fieldnames = ['aIDx', 'eIdx', 'tcIdx', 'swIdx', 'denIdx', 'dsnowIdx']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Integer')

        ## Write Double fields
        fieldnames = ['InitDensityScaling', 'ThermoDeltaTScaling', 'outputFreq', 'aSnow', 'aIce',
                      'cldFrac', 't0wet', 't0dry', 'K', 'adThresh', 'teThresh']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Double')

        ## Write DoubleMat fields
            ## mattype = 2
        fieldnames = ['aValue', 'teValue', 'dulwrfValue', 'szaValue', 'cotValue', 'ccsnowValue', 'cciceValue']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'DoubleMat', mattype = 2, timeserieslength = md.mesh.numberofelements + 1, yts = md.constants.yts)

            ## mattype = 3
        fieldnames = ['Dzini', 'Dini', 'Reini', 'Gdnini', 'Gspini', 'ECini', 'Wini', 'Aini', 'Adiffini', 'Tini']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'DoubleMat', mattype = 3)

        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'Sizeini', format = 'IntMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')

        if self.ismappedforcing:
            execute.WriteData(fid, prefix, obj = self, fieldname = 'mappedforcingpoint', format ='IntMat', mattype = 2)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'mappedforcingelevation', format ='DoubleMat', mattype = 3)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'lapseTaValue', format ='Double')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'lapsedlwrfValue', format ='Double')

        ## Calculate dt from forcings
        ## NOTE: Taken from $ISSM_DIR/src/m/classes/SMBgemb.py
        if (np.any(self.P[-1] - self.Ta[-1] != 0) | np.any(self.V[-1] - self.Ta[-1] != 0) | np.any(self.dswrf[-1] - self.Ta[-1] != 0) | np.any(self.dlwrf[-1] - self.Ta[-1] != 0) | np.any(self.eAir[-1] - self.Ta[-1] != 0) | np.any(self.pAir[-1] - self.Ta[-1] != 0)):
            raise IOError('All GEMB forcings (Ta, P, V, dswrf, dlwrf, eAir, pAir) must have the same time steps in the final row!')

        if ((np.ndim(self.teValue)>1) & np.any(self.teValue[-1] - self.Ta[-1] != 0)):
            raise IOError('If GEMB forcing teValue is transient, it must have the same time steps as input Ta in the final row!')
        if ((np.ndim(self.dswdiffrf)>1) & np.any(self.dswdiffrf[-1] - self.Ta[-1] != 0)):
            raise IOError('If GEMB forcing dswdiffrf is transient, it must have the same time steps as input Ta in the final row!')
        if ((np.ndim(self.aValue)>1) & np.any(self.aValue[-1] - self.Ta[-1] != 0)):
            raise IOError('If GEMB forcing aValue is transient, it must have the same time steps as input Ta in the final row!')
        if ((np.ndim(self.dulwrfValue)>1) & np.any(self.dulwrfValue[-1] - self.Ta[-1] != 0)):
            raise IOError('If GEMB forcing dulwrfValue is transient, it must have the same time steps as input Ta in the final row!')
        if ((np.ndim(self.szaValue)>1) & np.any(self.szaValue[-1] - self.Ta[-1] != 0)):
            raise IOError('If GEMB forcing szaValue is transient, it must have the same time steps as input Ta in the final row!')
        if ((np.ndim(self.cotValue)>1) & np.any(self.cotValue[-1] - self.Ta[-1] != 0)):
            raise IOError('If GEMB forcing cotValue is transient, it must have the same time steps as input Ta in the final row!')
        if ((np.ndim(self.ccsnowValue)>1) & np.any(self.ccsnowValue[-1] - self.Ta[-1] != 0)):
            raise IOError('If GEMB forcing ccsnowValue is transient, it must have the same time steps as input Ta in the final row!')
        if ((np.ndim(self.cciceValue)>1) & np.any(self.cciceValue[-1] - self.Ta[-1] != 0)):
            raise IOError('If GEMB forcing cciceValue is transient, it must have the same time steps as input Ta in the final row!')

        time = self.Ta[-1]  # Assume all forcings are on the same time step
        dtime = np.diff(time, n=1, axis=0)
        dt = min(dtime)

        execute.WriteData(fid, prefix, name = 'md.smb.dt', data = dt, format = 'Double', scale = md.constants.yts)
        
        # Check if smb_dt goes evenly into transient core time step
        if (md.timestepping.time_step % dt >= 1e-10):
            raise IOError('smb_dt/dt = {}. The number of SMB time steps in one transient core time step has to be an an integer'.format(md.timestepping.time_step / dt))
        # Make sure that adaptive time step is off
        if md.timestepping.__class__.__name__ == 'timesteppingadaptive':
            raise IOError('GEMB cannot be run with adaptive timestepping.  Check class type of md.timestepping')

        execute.WriteData(fid, prefix, name = 'md.smb.requested_outputs', data = self.process_outputs(md), format = 'StringArray')

## ------------------------------------------------------
## smb.gradients
## ------------------------------------------------------
@class_registry.register_class
class gradients(class_registry.manage_state):
    """
    Gradient-based surface mass balance model for ISSM.

    This class implements a gradient-based SMB model where SMB varies linearly with
    elevation relative to a reference elevation and SMB. Different gradients can be
    specified for accumulation and ablation regimes.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    href : ndarray, default=np.nan
        Reference elevation from which deviation is used to calculate SMB adjustment [m].
    smbref : ndarray, default=np.nan
        Reference SMB from which deviation is calculated [m/yr ice equiv].
    b_pos : ndarray, default=np.nan
        Slope of elevation-SMB regression line for accumulation regime.
    b_neg : ndarray, default=np.nan
        Slope of elevation-SMB regression line for ablation regime.
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : list, default=['default']
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the gradient SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the gradient SMB parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Notes
    -----
    SMB is calculated as:
    SMB = smbref + gradient * (elevation - href)
    where gradient = b_pos for positive SMB or b_neg for negative SMB.

    Examples
    --------
    md.smb = pyissm.param.smb.gradients()
    md.smb.href = reference_elevation
    md.smb.smbref = reference_smb
    md.smb.b_pos = positive_gradient
    md.smb.b_neg = negative_gradient
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.href = np.nan
        self.smbref = np.nan
        self.b_pos = np.nan
        self.b_neg = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'issmbgradients', 'is smb gradients method activated (0 or 1, default is 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'href', 'reference elevation from which deviation is used to calculate SMB adjustment in smb gradients method'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'smbref', 'reference smb from which deviation is calculated in smb gradients method [m/yr ice equiv]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'b_pos', 'slope of hs - smb regression line for accumulation regime required if smb gradients is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'b_neg', 'slope of hs - smb regression line for ablation regime required if smb gradients is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.gradients Class'
        return s
    
    # Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self,
                        md = None,
                        return_default_outputs = False):
        """
        Process requested outputs, expanding 'default' to appropriate outputs.

        Parameters
        ----------
        md : ISSM model object, optional
            Model object containing mesh information.
        return_default_outputs : bool, default=False
            Whether to also return the list of default outputs.
            
        Returns
        -------
        outputs : list
            List of output strings with 'default' expanded to actual output names.
        default_outputs : list, optional
            Returned only if `return_default_outputs=True`.
        """

        outputs = []

        ## Set default_outputs
        default_outputs = ['SmbMassBalance']

        ## Loop through all requested outputs
        for item in self.requested_outputs:
            
            ## Process default outputs
            if item == 'default':
                    outputs.extend(default_outputs)

            ## Append other requested outputs (not defaults)
            else:
                outputs.append(item)

        if return_default_outputs:
            return outputs, default_outputs
        return outputs
        
    # Marshall method for saving the smb.gradients parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [smb.gradients] parameters to a binary file.

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
        execute.WriteData(fid, prefix, name = 'md.smb.model', data = 6, format = 'Integer')

        execute.WriteData(fid, prefix, obj = self, fieldname = 'href', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'smbref', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'b_pos', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'b_neg', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.smb.requested_outputs', data = self.process_outputs(md), format = 'StringArray')

## ------------------------------------------------------
## smb.gradientscomponents
## ------------------------------------------------------
@class_registry.register_class
class gradientscomponents(class_registry.manage_state):
    """
    Component-based gradient surface mass balance model for ISSM.

    This class implements a gradient-based SMB model where accumulation and runoff
    components vary separately with elevation. Each component has its own reference
    value, reference elevation, and elevation gradient.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    accuref : ndarray, default=np.nan
        Reference value of the accumulation [m ice eq/yr].
    accualti : ndarray, default=np.nan
        Altitude at which the accumulation is equal to the reference value [m].
    accugrad : ndarray, default=np.nan
        Gradient of the variation of the accumulation (0 for uniform accumulation) [m ice eq/yr/m].
    runoffref : ndarray, default=np.nan
        Reference value of the runoff [m w.e. y^-1].
    runoffalti : ndarray, default=np.nan
        Altitude at which the runoff is equal to the reference value [m].
    runoffgrad : ndarray, default=np.nan
        Gradient of the variation of the runoff (0 for uniform runoff) [m w.e. m^-1 y^-1].
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : list, default=['default']
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the gradient components SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the gradient components SMB parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Notes
    -----
    SMB components are calculated as:
    accumulation = accuref + accugrad * (elevation - accualti)
    runoff = runoffref + runoffgrad * (elevation - runoffalti)
    SMB = accumulation - runoff

    Examples
    --------
    md.smb = pyissm.param.smb.gradientscomponents()
    md.smb.accuref = reference_accumulation
    md.smb.runoffref = reference_runoff
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.accuref = np.nan
        self.accualti = np.nan
        self.accugrad = np.nan
        self.runoffref = np.nan
        self.runoffalti = np.nan
        self.runoffgrad = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = ['default']

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'issmbgradients', 'is smb gradients method activated (0 or 1, default is 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'accuref', ' reference value of the accumulation'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'accualti', ' Altitude at which the accumulation is equal to the reference value'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'accugrad', ' Gradient of the variation of the accumulation (0 for uniform accumulation)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'runoffref', ' reference value of the runoff m w.e. y-1 (temperature times ddf)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'runoffalti', ' Altitude at which the runoff is equal to the reference value'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'runoffgrad', ' Gradient of the variation of the runoff (0 for uniform runoff) m w.e. m-1 y-1 (lapse rate times ddf)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.gradientscomponents Class'
        return s
    
    # Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self,
                        md = None,
                        return_default_outputs = False):
        """
        Process requested outputs, expanding 'default' to appropriate outputs.

        Parameters
        ----------
        md : ISSM model object, optional
            Model object containing mesh information.
        return_default_outputs : bool, default=False
            Whether to also return the list of default outputs.
            
        Returns
        -------
        outputs : list
            List of output strings with 'default' expanded to actual output names.
        default_outputs : list, optional
            Returned only if `return_default_outputs=True`.
        """

        outputs = []

        ## Set default_outputs
        default_outputs = ['SmbMassBalance']

        ## Loop through all requested outputs
        for item in self.requested_outputs:
            
            ## Process default outputs
            if item == 'default':

                    ## Add to default_outputs when steps_per_step > 1
                    if self.steps_per_step > 1:
                        default_outputs.append('SmbMassBalanceSubstep')

                    outputs.extend(default_outputs)

            ## Append other requested outputs (not defaults)
            else:
                outputs.append(item)

        if return_default_outputs:
            return outputs, default_outputs
        return outputs
        
    # Marshall method for saving the smb.gradientscomponents parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [smb.gradientscomponents] parameters to a binary file.

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
        execute.WriteData(fid, prefix, name = 'md.smb.model', data = 11, format = 'Integer')

        ## Write DoubleMat fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'accuref', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts, scale = 1. / md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'accugrad', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts, scale = 1. / md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'runoffref', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts, scale = 1. / md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'runoffgrad', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts, scale = 1. / md.constants.yts)
        
        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'accualti', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'runoffalti', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.smb.requested_outputs', data = self.process_outputs(md), format = 'StringArray')

## ------------------------------------------------------
## smb.gradientsela
## ------------------------------------------------------
@class_registry.register_class
class gradientsela(class_registry.manage_state):
    """
    Equilibrium Line Altitude (ELA) gradient surface mass balance model for ISSM.

    This class implements an ELA-based SMB model where SMB varies linearly with 
    elevation relative to the equilibrium line altitude. Different gradients are
    applied above and below the ELA, with optional caps on maximum and minimum SMB rates.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    ela : ndarray, default=np.nan
        Equilibrium line altitude from which deviation is used to calculate SMB [m a.s.l.].
    b_pos : ndarray, default=np.nan
        Vertical SMB gradient (dB/dz) above ELA [m ice eq./yr/m].
    b_neg : ndarray, default=np.nan
        Vertical SMB gradient (dB/dz) below ELA [m ice eq./yr/m].
    b_max : float, default=9999
        Upper cap on SMB rate [m ice eq./yr]. Default: 9999 (no cap).
    b_min : float, default=-9999
        Lower cap on SMB rate [m ice eq./yr]. Default: -9999 (no cap).
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : list, default=['default']
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the ELA gradient SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the ELA gradient SMB parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Notes
    -----
    SMB is calculated as:
    - For elevation > ELA: SMB = b_pos * (elevation - ELA)
    - For elevation < ELA: SMB = b_neg * (elevation - ELA)
    SMB is then clamped between b_min and b_max if specified.

    Examples
    --------
    md.smb = pyissm.param.smb.gradientsela()
    md.smb.ela = equilibrium_line_altitude
    md.smb.b_pos = positive_gradient  # Above ELA
    md.smb.b_neg = negative_gradient  # Below ELA
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.ela = np.nan
        self.b_pos = np.nan
        self.b_neg = np.nan
        self.b_max = 9999
        self.b_min = -9999
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = ['default']

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

        s += '\n   SMB gradients ela parameters:'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ela', ' equilibrium line altitude from which deviation is used to calculate smb using the smb gradients ela method [m a.s.l.]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'b_pos', ' vertical smb gradient (dB/dz) above ela'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'b_neg', ' vertical smb gradient (dB/dz) below ela'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'b_max', ' upper cap on smb rate, default: 9999 (no cap) [m ice eq./yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'b_min', ' lower cap on smb rate, default: -9999 (no cap) [m ice eq./yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.gradientsela Class'
        return s
    
    # Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self,
                        md = None,
                        return_default_outputs = False):
        """
        Process requested outputs, expanding 'default' to appropriate outputs.

        Parameters
        ----------
        md : ISSM model object, optional
            Model object containing mesh information.
        return_default_outputs : bool, default=False
            Whether to also return the list of default outputs.
            
        Returns
        -------
        outputs : list
            List of output strings with 'default' expanded to actual output names.
        default_outputs : list, optional
            Returned only if `return_default_outputs=True`.
        """

        outputs = []

        ## Set default_outputs
        default_outputs = ['SmbMassBalance']

        ## Loop through all requested outputs
        for item in self.requested_outputs:
            
            ## Process default outputs
            if item == 'default':
                    outputs.extend(default_outputs)

            ## Append other requested outputs (not defaults)
            else:
                outputs.append(item)

        if return_default_outputs:
            return outputs, default_outputs
        return outputs
        
    # Marshall method for saving the smb.gradientsela parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [smb.gradientsela] parameters to a binary file.

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
        execute.WriteData(fid, prefix, name = 'md.smb.model', data = 9, format = 'Integer')

        ## Write DoubleMat fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'ela', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'b_pos', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'b_neg', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'b_max', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'b_min', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        
        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.smb.requested_outputs', data = self.process_outputs(md), format = 'StringArray')

## ------------------------------------------------------
## smb.henning
## ------------------------------------------------------
@class_registry.register_class
class henning(class_registry.manage_state):
    """
    Henning surface mass balance model for ISSM.

    This class implements the Henning SMB parametrization, which is a specialized
    approach for modeling surface mass balance in ice sheet applications.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    smbref : ndarray, default=np.nan
        Reference surface mass balance [m/yr ice eq].
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : list, default=['default']
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the Henning SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the Henning SMB parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.smb = pyissm.param.smb.henning()
    md.smb.smbref = reference_smb_data
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.smbref = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = ['default']

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'mass_balance', 'surface mass balance [m/yr ice eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.henning Class'
        return s

    # Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self,
                        md = None,
                        return_default_outputs = False):
        """
        Process requested outputs, expanding 'default' to appropriate outputs.

        Parameters
        ----------
        md : ISSM model object, optional
            Model object containing mesh information.
        return_default_outputs : bool, default=False
            Whether to also return the list of default outputs.
            
        Returns
        -------
        outputs : list
            List of output strings with 'default' expanded to actual output names.
        default_outputs : list, optional
            Returned only if `return_default_outputs=True`.
        """

        outputs = []

        ## Set default_outputs
        default_outputs = ['SmbMassBalance']

        ## Loop through all requested outputs
        for item in self.requested_outputs:
            
            ## Process default outputs
            if item == 'default':
                    outputs.extend(default_outputs)

            ## Append other requested outputs (not defaults)
            else:
                outputs.append(item)

        if return_default_outputs:
            return outputs, default_outputs
        return outputs
        
    # Marshall method for saving the smb.henning parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [smb.henning] parameters to a binary file.

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
        execute.WriteData(fid, prefix, name = 'md.smb.model', data = 1, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'mass_balance', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.smb.requested_outputs', data = self.process_outputs(md), format = 'StringArray')

## ------------------------------------------------------
## smb.meltcomponents
## ------------------------------------------------------
@class_registry.register_class
class meltcomponents(class_registry.manage_state):
    """
    Melt component-based surface mass balance model for ISSM.

    This class implements a component-based SMB model that explicitly separates
    melt and refreeze processes. The surface mass balance is calculated as 
    SMB = accumulation - evaporation - melt + refreeze.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    accumulation : ndarray, default=np.nan
        Accumulated snow [m/yr ice eq].
    evaporation : ndarray, default=np.nan
        Amount of ice lost to evaporative processes [m/yr ice eq].
    melt : ndarray, default=np.nan
        Amount of ice melt in the ice column [m/yr ice eq].
    refreeze : ndarray, default=np.nan
        Amount of ice melt refrozen in the ice column [m/yr ice eq].
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : list, default=['default']
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the melt components SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the melt components SMB parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Notes
    -----
    The surface mass balance is computed as:
    SMB = accumulation - evaporation - melt + refreeze

    This formulation explicitly accounts for refreezing processes that can occur
    in firn layers, which is important for accurate SMB modeling in cold regions.

    Examples
    --------
    md.smb = pyissm.param.smb.meltcomponents()
    md.smb.accumulation = accumulation_data
    md.smb.melt = melt_data
    md.smb.refreeze = refreeze_data
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.accumulation = np.nan
        self.evaporation = np.nan
        self.melt = np.nan
        self.refreeze = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = ['default']

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters with melt (SMB = accumulation-evaporation-melt+refreeze):\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'accumulation', 'accumulated snow [m/yr ice eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'evaporation', 'mount of ice lost to evaporative processes [m/yr ice eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'melt', 'amount of ice melt in the ice column [m/yr ice eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'refreeze', 'amount of ice melt refrozen in the ice column [m/yr ice eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.meltcomponents Class'
        return s

    # Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self,
                        md = None,
                        return_default_outputs = False):
        """
        Process requested outputs, expanding 'default' to appropriate outputs.

        Parameters
        ----------
        md : ISSM model object, optional
            Model object containing mesh information.
        return_default_outputs : bool, default=False
            Whether to also return the list of default outputs.
            
        Returns
        -------
        outputs : list
            List of output strings with 'default' expanded to actual output names.
        default_outputs : list, optional
            Returned only if `return_default_outputs=True`.
        """

        outputs = []

        ## Set default_outputs
        default_outputs = ['SmbMassBalance']

        ## Loop through all requested outputs
        for item in self.requested_outputs:
            
            ## Process default outputs
            if item == 'default':
                    outputs.extend(default_outputs)

            ## Append other requested outputs (not defaults)
            else:
                outputs.append(item)

        if return_default_outputs:
            return outputs, default_outputs
        return outputs
        
    # Marshall method for saving the smb.meltcomponents parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [smb.meltcomponents] parameters to a binary file.

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
        execute.WriteData(fid, prefix, name = 'md.smb.model', data = 3, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'accumulation', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'evaporation', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'melt', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'refreeze', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.smb.requested_outputs', data = self.process_outputs(md), format = 'StringArray')

## ------------------------------------------------------
## smb.pdd
## ------------------------------------------------------
@class_registry.register_class
class pdd(class_registry.manage_state):
    """
    Positive Degree Day surface mass balance model for ISSM.

    This class implements a positive degree day (PDD) SMB model that calculates
    surface mass balance based on temperature and precipitation data. It supports
    multiple temperature and precipitation data sources, including delta-18-O and
    MUNGSM parametrizations for paleoclimate applications.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    precipitation : ndarray, default=np.nan
        Monthly surface precipitation [m/yr water eq].
    monthlytemperatures : ndarray, default=np.nan
        Monthly surface temperatures [K].
    desfac : float, default=0.5
        Desertification elevation factor (between 0 and 1) [m].
    s0p : ndarray, default=np.nan
        Elevation from precipitation source (between 0 and a few 1000s m) [m].
    s0t : ndarray, default=np.nan
        Elevation from temperature source (between 0 and a few 1000s m) [m].
    rlaps : float, default=6.5
        Present day lapse rate [degree/km].
    rlapslgm : float, default=6.5
        LGM lapse rate [degree/km].
    Pfac : ndarray, default=np.nan
        Time interpolation parameter for precipitation, 1D(year).
    Tdiff : ndarray, default=np.nan
        Time interpolation parameter for temperature, 1D(year).
    sealev : ndarray, default=np.nan
        Sea level [m], 1D(year).
    isdelta18o : int, default=0
        Is temperature and precipitation delta18o parametrisation activated (0 or 1).
    ismungsm : int, default=0
        Is temperature and precipitation mungsm parametrisation activated (0 or 1).
    issetpddfac : int, default=0
        Is user passing in defined PDD factors (0 or 1).
    delta18o : float, default=0
        Delta-18-O values [per mil].
    delta18o_surface : ndarray, default=np.nan
        Surface elevation of the delta18o site [m].
    temperatures_presentday : ndarray, default=np.nan
        Monthly present day surface temperatures [K].
    temperatures_lgm : ndarray, default=np.nan
        Monthly LGM surface temperatures [K].
    precipitations_presentday : ndarray, default=np.nan
        Monthly present day surface precipitation [m/yr water eq].
    precipitations_lgm : ndarray, default=np.nan
        Monthly LGM surface precipitation [m/yr water eq].
    pddfac_snow : ndarray, default=np.nan
        PDD factor for snow [mm ice equiv/day/degree C].
    pddfac_ice : ndarray, default=np.nan
        PDD factor for ice [mm ice equiv/day/degree C].
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : list, default=['default']
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the PDD SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the PDD SMB parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Notes
    -----
    The PDD model calculates melt based on the number of positive degree days,
    which is the sum of temperatures above freezing over a given time period.
    This approach is widely used in glaciology for its simplicity and effectiveness.

    Examples
    --------
    md.smb = pyissm.param.smb.pdd()
    md.smb.monthlytemperatures = temperature_data
    md.smb.precipitation = precipitation_data
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.precipitation = np.nan
        self.monthlytemperatures = np.nan
        self.desfac = 0.5
        self.s0p = np.nan
        self.s0t = np.nan
        self.rlaps = 6.5
        self.rlapslgm = 6.5
        self.Pfac = np.nan
        self.Tdiff = np.nan
        self.sealev = np.nan
        self.isdelta18o = 0
        self.ismungsm = 0
        self.issetpddfac = 0
        self.delta18o = 0
        self.delta18o_surface = np.nan
        self.temperatures_presentday = np.nan
        self.temperatures_lgm = np.nan
        self.precipitations_presentday = np.nan
        self.precipitations_lgm = np.nan
        self.pddfac_snow = np.nan
        self.pddfac_ice = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'isdelta18o', 'is temperature and precipitation delta18o parametrisation activated (0 or 1, default is 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ismungsm', 'is temperature and precipitation mungsm parametrisation activated (0 or 1, default is 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'issetpddfac', 'is user passing in defined pdd factors (0 or 1, default is 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'desfac', 'desertification elevation factor (between 0 and 1, default is 0.5) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 's0p', 'should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 's0t', 'should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rlaps', 'present day lapse rate [degree/km]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rlapslgm', 'LGM lapse rate [degree/km]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'monthlytemperatures', 'monthly surface temperatures [K], required if pdd is activated and delta18o not activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'precipitation', 'monthly surface precipitation [m/yr water eq], required if pdd is activated and delta18o or mungsm not activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'delta18o', 'delta18o [per mil], required if pdd is activated and delta18o activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'delta18o_surface', 'surface elevation of the delta18o site, required if pdd is activated and delta18o activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'temperatures_presentday', 'monthly present day surface temperatures [K], required if delta18o/mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'temperatures_lgm', 'monthly LGM surface temperatures [K], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'precipitations_presentday', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'precipitations_lgm', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Tdiff', 'time interpolation parameter for temperature, 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sealev', 'sea level [m], 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'temperatures_presentday', 'monthly present day surface temperatures [K], required if delta18o/mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'temperatures_lgm', 'monthly LGM surface temperatures [K], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'precipitations_presentday', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'precipitations_lgm', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Pfac', 'time interpolation parameter for precipitation, 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Tdiff', 'time interpolation parameter for temperature, 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sealev', 'sea level [m], 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.pdd Class'
        return s

    # Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self,
                        md = None,
                        return_default_outputs = False):
        """
        Process requested outputs, expanding 'default' to appropriate outputs.

        Parameters
        ----------
        md : ISSM model object, optional
            Model object containing mesh information.
        return_default_outputs : bool, default=False
            Whether to also return the list of default outputs.
            
        Returns
        -------
        outputs : list
            List of output strings with 'default' expanded to actual output names.
        default_outputs : list, optional
            Returned only if `return_default_outputs=True`.
        """

        outputs = []

        ## Set default_outputs
        default_outputs = ['SmbMassBalance']

        ## Loop through all requested outputs
        for item in self.requested_outputs:
            
            ## Process default outputs
            if item == 'default':
                    outputs.extend(default_outputs)

            ## Append other requested outputs (not defaults)
            else:
                outputs.append(item)

        if return_default_outputs:
            return outputs, default_outputs
        return outputs
        
    # Marshall method for saving the smb.pdd parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [smb.pdd] parameters to a binary file.

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
        execute.WriteData(fid, prefix, name = 'md.smb.model', data = 4, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'isdelta18o', format = 'Boolean')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'ismungsm', format = 'Boolean')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'issetpddfac', format = 'Boolean')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'desfac', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 's0p', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 's0t', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rlaps', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rlapslgm', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.smb.requested_outputs', data = self.process_outputs(md), format = 'StringArray')

        ## Write conditional fields
        if (self.isdelta18o == 0 and self.ismungsm == 0):
            execute.WriteData(fid, prefix, obj = self, fieldname = 'monthlytemperatures', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'precipitation', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        elif self.isdelta18o:
            execute.WriteData(fid, prefix, obj = self, fieldname = 'temperatures_presentday', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'temperatures_lgm', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'precipitations_presentday', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'precipitations_lgm', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'delta18o_surface', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'delta18o', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'Tdiff', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'sealev', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)
        elif self.ismungsm:
            execute.WriteData(fid, prefix, obj = self, fieldname = 'temperatures_presentday', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'temperatures_lgm', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'precipitations_presentday', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'precipitations_lgm', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'Pfac', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'Tdiff', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'sealev', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)

        if self.issetpddfac:
            execute.WriteData(fid, prefix, ovj = self, fieldname = 'pddfac_snow', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute.WriteData(fid, prefix, ovj = self, fieldname = 'pddfac_ice', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

## ------------------------------------------------------
## smb.pddSicopolis
## ------------------------------------------------------
@class_registry.register_class
class pddSicopolis(class_registry.manage_state):
    """
    SICOPOLIS-style Positive Degree Day surface mass balance model for ISSM.

    This class implements the SICOPOLIS PDD scheme (Calov & Greve, 2005) for 
    surface mass balance calculations. It includes temperature and precipitation
    anomalies, firn warming effects, and desertification corrections based on
    the SICOPOLIS ice sheet model approach.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values 
        in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    precipitation : ndarray, default=np.nan
        Monthly surface precipitation [m/yr water eq].
    monthlytemperatures : ndarray, default=np.nan
        Monthly surface temperatures [K].
    temperature_anomaly : ndarray, default=np.nan
        Anomaly to monthly reference temperature (additive) [K].
    precipitation_anomaly : ndarray, default=np.nan
        Anomaly to monthly precipitation (multiplicative, e.g. q = q0*exp(0.070458*DeltaT)) [unitless].
    smb_corr : ndarray, default=np.nan
        Correction of SMB after PDD call [m/a].
    desfac : float, default=-np.log(2.0)/1000
        Desertification elevation factor. Default: -log(2.0)/1000.
    s0p : ndarray, default=np.nan
        Elevation from precipitation source (between 0 and a few 1000s m) [m].
    s0t : ndarray, default=np.nan
        Elevation from temperature source (between 0 and a few 1000s m) [m].
    rlaps : float, default=7.4
        Present day lapse rate [degree/km]. Default: 7.4.
    isfirnwarming : int, default=1
        Is firn warming (Reeh 1991) activated (0 or 1). Default: 1.
    steps_per_step : int, default=1
        Number of SMB steps per time step.
    averaging : int, default=0
        Averaging method from short to long steps. 0: Arithmetic, 1: Geometric, 2: Harmonic.
    requested_outputs : list, default=['default']
        Additional outputs requested (TemperaturePDD, SmbAccumulation, SmbMelt).

    Methods
    -------
    __init__(self, other=None)
        Initializes the SICOPOLIS PDD SMB parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the SICOPOLIS PDD SMB parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Notes
    -----
    This implementation follows the SICOPOLIS PDD scheme as described in:
    Calov, R., & Greve, R. (2005). A semi-analytical solution for the positive 
    degree-day model with stochastic temperature variations. Journal of Glaciology, 
    51(172), 173-175.

    The firn warming correction (Reeh, 1991) adjusts melt rates based on firn
    temperature, which is important for accurate SMB calculations in accumulation zones.

    Examples
    --------
    md.smb = pyissm.param.smb.pddSicopolis()
    md.smb.monthlytemperatures = temperature_data
    md.smb.precipitation = precipitation_data
    md.smb.temperature_anomaly = temp_anomaly
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.precipitation = np.nan
        self.monthlytemperatures = np.nan
        self.temperature_anomaly = np.nan
        self.precipitation_anomaly = np.nan
        self.smb_corr = np.nan
        self.desfac = -np.log(2.0) / 1000
        self.s0p = np.nan
        self.s0t = np.nan
        self.rlaps = 7.4
        self.isfirnwarming = 1
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = ['default']

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'
        s += '   SICOPOLIS PDD scheme (Calov & Greve, 2005):\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'monthlytemperatures', 'monthly surface temperatures [K]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'precipitation', 'monthly surface precipitation [m/yr water eq]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'temperature_anomaly', 'anomaly to monthly reference temperature (additive [K])'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'precipitation_anomaly', 'anomaly to monthly precipitation (multiplicative, e.g. q = q0*exp(0.070458*DeltaT) after Huybrechts (2002)) [no unit])'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'smb_corr', 'correction of smb after PDD call [m/a]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 's0p', 'should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 's0t', 'should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rlaps', 'present day lapse rate (default is 7.4 degree/km)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'desfac', 'desertification elevation factor (default is -log(2.0)/1000)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isfirnwarming', 'is firnwarming (Reeh 1991) activated (0 or 1, default is 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested (TemperaturePDD, SmbAccumulation, SmbMelt)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.pddSicopolis Class'
        return s

    # Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self,
                        md = None,
                        return_default_outputs = False):
        """
        Process requested outputs, expanding 'default' to appropriate outputs.

        Parameters
        ----------
        md : ISSM model object, optional
            Model object containing mesh information.
        return_default_outputs : bool, default=False
            Whether to also return the list of default outputs.
            
        Returns
        -------
        outputs : list
            List of output strings with 'default' expanded to actual output names.
        default_outputs : list, optional
            Returned only if `return_default_outputs=True`.
        """

        outputs = []

        ## Set default_outputs
        default_outputs = ['SmbMassBalance']

        ## Loop through all requested outputs
        for item in self.requested_outputs:
            
            ## Process default outputs
            if item == 'default':
                    outputs.extend(default_outputs)

            ## Append other requested outputs (not defaults)
            else:
                outputs.append(item)

        if return_default_outputs:
            return outputs, default_outputs
        return outputs
        
    # Marshall method for saving the smb.pddSicopolis parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [smb.pddSicopolis] parameters to a binary file.

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
        execute.WriteData(fid, prefix, name = 'md.smb.model', data = 10, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'isfirnwarming', format = 'Boolean')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'desfac', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 's0p', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 's0t', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'rlaps', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'monthlytemperatures', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'precipitation', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'temperature_anomaly', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'precipitation_anomaly', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'smb_corr', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.smb.requested_outputs', data = self.process_outputs(md), format = 'StringArray')
