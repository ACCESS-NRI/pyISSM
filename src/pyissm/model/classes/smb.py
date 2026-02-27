import numpy as np
import warnings
from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm.model import execute, mesh
from pyissm import model

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
    _process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.smb = pyissm.model.classes.smb.default()
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

        s += '{}\n'.format(class_utils._field_display(self, 'mass_balance', 'surface mass balance [m/yr ice eq]'))
        s += '{}\n'.format(class_utils._field_display(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(class_utils._field_display(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(class_utils._field_display(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.default Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [smb.default] fields to 3D
        """
        self.mass_balance = mesh._project_3d(md, vector = self.mass_balance, type = 'node')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        # Early return if required analysis/solutions are not present
        if solution == 'TransientSolution' and not md.transient.issmb:
            return
        
        if 'MasstransportAnalysis' in analyses:
            class_utils._check_field(md, fieldname = 'smb.mass_balance', timeseries= True, allow_nan = False, allow_inf = False)
        if 'BalancethicknessAnalysis' in analyses:
            class_utils._check_field(md, fieldname = 'smb.mass_balance', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.steps_per_step', ge = 1, scalar = True)
        class_utils._check_field(md, fieldname = 'smb.averaging', scalar = True, values = [0, 1, 2])
        class_utils._check_field(md, fieldname = 'smb.requested_outputs', string_list = True)
        
        return md
    
    # Initialise empty fields of correct dimensions
    def _initialize(self, md):
        """
        Initialise empty fields in smb.default.
        """

        if np.all(np.isnan(self.mass_balance)):
            self.mass_balance = np.zeros((md.mesh.numberofvertices))
            warnings.warn('pyissm.model.classes.smb.default: smb.mass_balance not specified -- set to 0.')

        return self
    
    # Process requested outputs, expanding 'default' to appropriate outputs
    def _process_outputs(self,
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
        execute._write_model_field(fid, prefix, name = 'md.smb.model', data = 1, format = 'Integer')

        ## Write fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'mass_balance', format = 'DoubleMat', scale = 1. / md.constants.yts, mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, name = 'md.smb.requested_outputs', data = self._process_outputs(md), format = 'StringArray')

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
    _process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.smb = pyissm.model.classes.smb.arma()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.num_basins = 0
        self.num_params = 0
        self.num_breaks = 0
        self.arma_timestep = 0
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

        s += '{}\n'.format(class_utils._field_display(self, 'num_basins', 'number of different basins [unitless]'))
        s += '{}\n'.format(class_utils._field_display(self, 'basin_id', 'basin number assigned to each element [unitless]'))
        s += '{}\n'.format(class_utils._field_display(self, 'num_breaks', 'number of different breakpoints in the piecewise-polynomial (separating num_breaks+1 periods)'))
        s += '{}\n'.format(class_utils._field_display(self, 'num_params', 'number of different parameters in the piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.)'))
        s += '{}\n'.format(class_utils._field_display(self, 'polynomialparams', 'coefficients for the polynomial (const,trend,quadratic,etc.),dim1 for basins,dim2 for periods,dim3 for orders, ex: polyparams=cat(num_params,intercepts,trendlinearcoefs,trendquadraticcoefs)'))
        s += '{}\n'.format(class_utils._field_display(self, 'datebreaks', 'dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr]'))
        s += '{}\n'.format(class_utils._field_display(self, 'ar_order', 'order of the autoregressive model [unitless]'))
        s += '{}\n'.format(class_utils._field_display(self, 'ma_order', 'order of the moving-average model [unitless]'))
        s += '{}\n'.format(class_utils._field_display(self, 'arma_timestep', 'time resolution of the ARMA model [yr]'))
        s += '{}\n'.format(class_utils._field_display(self, 'arlag_coefs', 'basin-specific vectors of AR lag coefficients [unitless]'))
        s += '{}\n'.format(class_utils._field_display(self, 'malag_coefs', 'basin-specific vectors of MA lag coefficients [unitless]'))
        s += '{}\n'.format(class_utils._field_display(self, 'lapserates', 'basin-specific SMB lapse rates applied in each elevation bin, 1 row per basin, 1 column per bin, dimension 3 can be of size 12 to prescribe monthly varying values [m ice eq yr^-1 m^-1] (default: no lapse rate)'))
        s += '{}\n'.format(class_utils._field_display(self, 'elevationbins', 'basin-specific separations between elevation bins, 1 row per basin, 1 column per limit between bins, dimension 3 can be of size 12 to prescribe monthly varying values [m] (default: no basin separation)'))
        s += '{}\n'.format(class_utils._field_display(self, 'refelevation', 'basin-specific reference elevations at which SMB is calculated, and from which SMB is downscaled using lapserates (default: basin mean elevation) [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(class_utils._field_display(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(class_utils._field_display(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.arma Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [smb.arma]fields to 3D
        """
        warnings.warn('pyissm.model.classes.smb.arma._extrude: 3D extrusion not implemented for smb.arma. Returning unchanged (2D) smb fields.')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        if 'MasstransportAnalysis' in analyses:
            nbas = md.smb.num_basins
            nprm = md.smb.num_params
            nbrk = md.smb.num_breaks
            class_utils._check_field(md, fieldname = 'smb.num_basins', scalar = True, gt = 0, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.num_params', scalar = True, gt = 0, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.num_breaks', scalar = True, ge = 0, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.basin_id', ge = 0, le = md.smb.num_basins, size = (md.mesh.numberofelements, ), allow_inf = False)

            if nbas > 1 and nbrk >= 1 and nprm > 1:
                class_utils._check_field(md, fieldname = 'smb.polynomialparams', size = (nbas, nbrk + 1, nprm), numel = nbas * (nbrk + 1) * nprm, allow_nan = False, allow_inf = False)
            elif nbas == 1:
                class_utils._check_field(md, fieldname = 'smb.polynomialparams', size = (nprm, nbrk + 1), numel = nbas * (nbrk + 1) * nprm, allow_nan = False, allow_inf = False)
            elif nbrk == 0:
                class_utils._check_field(md, fieldname = 'smb.polynomialparams', size = (nbas, nprm), numel = nbas * (nbrk + 1) * nprm, allow_nan = False, allow_inf = False)
            elif nprm == 1:
                class_utils._check_field(md, fieldname = 'smb.polynomialparams', size = (nbas, nbrk), numel = nbas * (nbrk + 1) * nprm, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.ar_order', scalar = True, ge = 0, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.ma_order', scalar = True, ge = 0, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.arma_timestep', scalar = True, ge = md.timestepping.time_step, allow_nan = False, allow_inf = False) # Autoregression time step cannot be finer than ISSM timestep
            class_utils._check_field(md, fieldname = 'smb.arlag_coefs', size = (nbas, md.smb.ar_order), allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.malag_coefs', size = (nbas, md.smb.ma_order), allow_nan = False, allow_inf = False)
            if nbrk > 0:
                class_utils._check_field(md, fieldname = 'smb.datebreaks', size = (nbas, nbrk), allow_nan = False, allow_inf = False)
            elif np.size(md.smb.datebreaks) == 0 or np.all(np.isnan(md.smb.datebreaks)):
                pass
            else:
                raise RuntimeError('pyissm.model.classes.smb.arma.check_consistency: md.smb.num_breaks is 0 but md.smb.datebreaks is not empty')

            if np.any(np.isnan(self.refelevation) is False) or np.size(self.refelevation) > 1:
                if len(np.shape(self.refelevation)) == 1:
                    self.refelevation = np.array([self.refelevation])
                class_utils._check_field(md, fieldname = 'smb.refelevation', size = (1, nbas), numel = nbas, ge = 0, allow_nan = False, allow_inf = False)

            if (np.any(np.isnan(self.lapserates) is False) or np.size(self.lapserates) > 1):
                nbas = md.smb.num_basins
                if len(np.shape(self.lapserates)) == 1:
                    nbins = 1
                    self.lapserates = np.reshape(self.lapserates,[nbas,nbins,1])
                elif(len(np.shape(self.lapserates)) == 2):
                    nbins = np.shape(self.lapserates)[1]
                    self.lapserates = np.reshape(self.lapserates,[nbas,nbins,1])
                elif(len(np.shape(self.lapserates)) == 3):
                    nbins = np.shape(self.lapserates)[1]
                ntmlapse = np.shape(self.lapserates)[2]
                if len(np.shape(self.elevationbins)) < 3:
                    self.elevationbins = np.reshape(self.elevationbins,[nbas,max(1,nbins-1),ntmlapse])
                class_utils._check_field(md, fieldname = 'smb.lapserates', size = (nbas, nbins, ntmlapse), numel = md.smb.num_basins * nbins * ntmlapse, allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.elevationbins', size = (nbas,max(1,nbins-1),ntmlapse), numel = nbas*max(1,nbins-1)*ntmlapse, allow_nan = False, allow_inf = False)
                for rr in range(nbas):
                    if(np.all(self.elevationbins[rr,0:-1]<=self.elevationbins[rr,1:])==False):
                        raise TypeError('pyissm.model.classes.smb.arma.check_consistency: md.smb.elevationbins should have rows in order of increasing elevation')
            elif (np.any(np.isnan(self.elevationbins) is False) or np.size(self.elevationbins) > 1):
                # Elevationbins specified but not lapserates: this will inevitably lead to inconsistencies
                nbas = md.smb.num_basins
                if len(np.shape(self.elevationbins)) == 1:
                    nbins = 1
                    self.elevationbins = np.reshape(self.elevationbins,[nbas,nbins,1])
                elif(len(np.shape(self.lapserates)) == 2):
                    nbins = np.shape(self.elevationbins)[1]
                    self.elevationbins = np.reshape(self.elevationbins,[nbas,nbins,1])
                elif(len(np.shape(self.lapserates)) == 3):
                    nbins = np.shape(self.lapserates)[1]
                nbins = nbins - 1
                ntmlapse = np.shape(self.lapserates)[2]
                class_utils._check_field(md, fieldname = 'smb.lapserates', size = (nbas, nbins * ntmlapse), numel = nbas * nbins * ntmlapse, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.elevationbins', size = (nbas, max(1, nbins - 1) * ntmlapse), numel = nbas * max(1, nbins - 1) * ntmlapse, allow_nan = False, allow_inf = False)

        class_utils._check_field(md, fieldname = 'smb.steps_per_step', ge = 1, scalar = True)
        class_utils._check_field(md, fieldname = 'smb.averaging', scalar = True, values = [0, 1, 2])
        class_utils._check_field(md, fieldname = 'smb.requested_outputs', string_list = True)

        return md
    
    # Initialise empty fields of correct dimensions
    def _initialize(self, md):
        """
        Initialise empty fields in smb.arma.
        """
                
        if self.ar_order == 0:
            self.ar_order = 1 # Dummy 1 value for autoregression
            self.arlag_coefs = np.zeros((self.num_basins, self.ar_order)) # Autoregression coefficients all set to 0
            warnings.warn('pyissm.model.classes.smb.arma: smb.ar_order (order of autoregressive model) not specified -- order of autoregressive model set to 0.')

        if self.ma_order == 0:
            self.ma_order = 1 # Dummy 1 value for moving-average
            self.malag_coefs = np.zeros((self.num_basins, self.ma_order)) # Moving-average coefficients all set to 0
            warnings.warn('pyissm.model.classes.smb.arma: smb.ma_order (order of moving-average model) not specified -- order of moving-average model set to 0.')

        if self.arma_timestep == 0:
            self.arma_timestep = md.timestepping.time_step # ARMA model has no prescribed time step
            warnings.warn('pyissm.model.classes.smb.arma: smb.arma_timestep (timestep of ARMA model) not specified -- set to md.timestepping.time_step.')

        if np.all(np.isnan(self.arlag_coefs)):
            self.arlag_coefs = np.zeros((self.num_basins, self.ar_order)) # Autoregression model of order 0
            warnings.warn('pyissm.model.classes.smb.arma: smb.arlag_coefs (AR lag coefficients) not specified -- order of autoregressive model set to 0.')

        if np.all(np.isnan(self.malag_coefs)):
            self.malag_coefs = np.zeros((self.num_basins, self.ma_order)) # Moving-average model of order 0
            warnings.warn('pyissm.model.classes.smb.arma: smb.malag_coefs (MA lag coefficients) not specified -- order of moving-average model set to 0.')

        return self

    # Process requested outputs, expanding 'default' to appropriate outputs
    def _process_outputs(self,
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
        nper = self.num_breaks + 1
        polyParams_scaled_2d = np.zeros((self.num_basins, nper * self.num_params))
        if self.num_params > 1:
            # Case 3D
            if self.num_basins > 1 and nper > 1:
                for ii in range(self.num_params):
                    polyParams_scaled[:, :, ii] = polyParams_scaled[:, :, ii] * (1 / md.constants.yts) ** (ii + 1)
                # Fit in 2D array
                for ii in range(self.num_params):
                    polyParams_scaled_2d[:, ii * nper : (ii + 1) * nper] = 1 * polyParams_scaled[:, :, ii]
            # Case 2D and higher-order params at increasing row index
            elif self.num_basins == 1:
                for ii in range(self.num_params):
                    polyParams_scaled[ii, :] = polyParams_scaled[ii, :] * (1 / md.constants.yts) ** (ii + 1)
                # Fit in row array
                for ii in range(self.num_params):
                    polyParams_scaled_2d[0, ii * nper : (ii + 1) * nper] = 1 * polyParams_scaled[ii, :]
            # Case 2D and higher-order params at increasing column index
            elif nper == 1:
                for ii in range(self.num_params):
                    polyParams_scaled[:, ii] = polyParams_scaled[:, ii] * (1 / md.constants.yts) ** (ii + 1)
                # 2D array is already in correct format
                polyParams_scaled_2d = np.copy(polyParams_scaled)
        else:
            polyParams_scaled   = polyParams_scaled * (1 / md.constants.yts)
            # 2D array is already in correct format
            polyParams_scaled_2d = np.copy(polyParams_scaled)

        if nper == 1:
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
        execute._write_model_field(fid, prefix, name = 'md.smb.model', data = 13, format = 'Integer')

        ## Write Integer fields
        fieldnames = ['num_basins', 'num_breaks', 'num_params', 'ar_order', 'ma_order', 'steps_per_step', 'averaging']
        for field in fieldnames:
            execute._write_model_field(fid, prefix, obj = self, fieldname = field, format = 'Integer')
        execute._write_model_field(fid, prefix, name = 'md.smb.num_bins', data = nbins, format = 'Integer')

        ## Write DoubleMat fields
        execute._write_model_field(fid, prefix, name = 'md.smb.polynomialparams', data = polyParams_scaled_2d,  format = 'DoubleMat')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'arlag_coefs', format = 'DoubleMat', yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'malag_coefs', format = 'DoubleMat', yts = md.constants.yts)
        execute._write_model_field(fid, prefix, name = 'md.smb.datebreaks', data = dbreaks, format = 'DoubleMat', scale = md.constants.yts)
        execute._write_model_field(fid, prefix, name = 'md.smb.lapserates', data = temp_lapse_rates_2d, format = 'DoubleMat', scale = 1. / md.constants.yts, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, name = 'md.smb.elevationbins', data = temp_elevation_bins_2d, format = 'DoubleMat')
        execute._write_model_field(fid, prefix, name = 'md.smb.refelevation', data = temp_ref_elevation, format = 'DoubleMat')

        ## Write other fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'arma_timestep', format = 'Double', scale = md.constants.yts)
        execute._write_model_field(fid, prefix, name = 'md.smb.basin_id', data = self.basin_id - 1, format = 'IntMat', mattype = 2)  # 0-indexed
        execute._write_model_field(fid, prefix, name = 'md.smb.requested_outputs', data = self._process_outputs(md), format = 'StringArray')

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
    _process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Notes
    -----
    The surface mass balance is computed as:
    SMB = accumulation - runoff - evaporation

    Examples
    --------
    md.smb = pyissm.model.classes.smb.components()
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

        s += '{}\n'.format(class_utils._field_display(self, 'accumulation', 'accumulated snow [m/yr ice eq]'))
        s += '{}\n'.format(class_utils._field_display(self, 'runoff', 'amount of ice melt lost from the ice column [m/yr ice eq]'))
        s += '{}\n'.format(class_utils._field_display(self, 'evaporation', 'mount of ice lost to evaporative processes [m/yr ice eq]'))
        s += '{}\n'.format(class_utils._field_display(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(class_utils._field_display(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(class_utils._field_display(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.components Class'
        return s

    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [smb.components]fields to 3D
        """
        self.accumulation = mesh._project_3d(md, vector = self.accumulation, type = 'node')
        self.runoff = mesh._project_3d(md, vector = self.runoff, type = 'node')
        self.evaporation = mesh._project_3d(md, vector = self.evaporation, type = 'node')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        if 'MasstransportAnalysis' in analyses:
            class_utils._check_field(md, fieldname = 'smb.accumulation', timeseries = True, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.runoff', timeseries = True, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.evaporation', timeseries = True, allow_nan = False, allow_inf = False)
        if 'BalancethicknessAnalysis' in analyses:
            class_utils._check_field(md, fieldname = 'smb.accumulation', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.runoff', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.evaporation', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.steps_per_step', ge = 1, scalar = True)
        class_utils._check_field(md, fieldname = 'smb.averaging', scalar = True, values = [0, 1, 2])
        class_utils._check_field(md, fieldname = 'smb.requested_outputs', string_list = True)

        return md

    # Initialise empty fields of correct dimensions
    def _initialize(self, md):
        """
        Initialise empty fields in smb.components.
        """

        if np.all(np.isnan(self.accumulation)):
            self.accumulation = np.zeros((md.mesh.numberofvertices))
            warnings.warn('pyissm.model.classes.smb.components: no SMB.accumulation specified -- values set as 0.')

        if np.all(np.isnan(self.evaporation)):
            self.evaporation = np.zeros((md.mesh.numberofvertices))
            warnings.warn('pyissm.model.classes.smb.components: no SMB.evaporation specified -- values set as 0.')

        if np.all(np.isnan(self.runoff)):
            self.runoff = np.zeros((md.mesh.numberofvertices))
            warnings.warn('pyissm.model.classes.smb.components: no SMB.runoff specified -- values set as 0.')

        return self
    
    # Process requested outputs, expanding 'default' to appropriate outputs
    def _process_outputs(self,
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
        execute._write_model_field(fid, prefix, name = 'md.smb.model', data = 2, format = 'Integer')

        ## Write fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'accumulation', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'runoff', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'evaporation', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute._write_model_field(fid, prefix, name = 'md.smb.requested_outputs', data = self._process_outputs(md), format = 'StringArray')


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
    _process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.smb = pyissm.model.classes.smb.d18opdd()
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
        self.requested_outputs = ['default']

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

        s += '{}\n'.format(class_utils._field_display(self, 'isd18opd', 'is delta18o parametrisation from present day temperature and precipitation activated (0 or 1, default is 0)'))
        s += '{}\n'.format(class_utils._field_display(self, 'issetpddfac', 'is user passing in defined pdd factors (0 or 1, default is 0)'))
        s += '{}\n'.format(class_utils._field_display(self, 'desfac', 'desertification elevation factor (between 0 and 1, default is 0.5) [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 's0p', 'should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 's0t', 'should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'rlaps', 'present day lapse rate [degree/km]'))
        s += '{}\n'.format(class_utils._field_display(self, 'temperatures_presentday', 'monthly present day surface temperatures [K], required if delta18o/mungsm is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'precipitations_presentday', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'istemperaturescaled', 'if delta18o parametrisation from present day temperature and precipitation is activated, is temperature scaled to delta18o value? (0 or 1, default is 1)'))
        s += '{}\n'.format(class_utils._field_display(self, 'isprecipscaled', 'if delta18o parametrisation from present day temperature and precipitation is activated, is precipitation scaled to delta18o value? (0 or 1, default is 1)'))
        s += '{}\n'.format(class_utils._field_display(self, 'temperatures_reconstructed', 'monthly historical surface temperatures [K], required if delta18o/mungsm/d18opd is activated and istemperaturescaled is not activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'precipitations_reconstructed', 'monthly historical precipitation [m/yr water eq], required if delta18o/mungsm/d18opd is activated and isprecipscaled is not activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'delta18o', 'delta18o [per mil], required if pdd is activated and delta18o activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'dpermil', 'degree per mil, required if d18opd is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'f', 'precip/temperature scaling factor, required if d18opd is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'pddfac_snow', 'Pdd factor for snow for all the domain [mm ice equiv/day/degree C]'))
        s += '{}\n'.format(class_utils._field_display(self, 'pddfac_ice', 'Pdd factor for ice for all the domain [mm ice equiv/day/degree C]'))
        s += '{}\n'.format(class_utils._field_display(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(class_utils._field_display(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(class_utils._field_display(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.d18opdd Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [smb.d18opdd] fields to 3D
        """
        if self.isd18opd:
            self.temperatures_presentday = mesh._project_3d(md, vector = self.temperatures_presentday, type = 'node')
        if self.isd18opd:
            self.precipitations_presentday = mesh._project_3d(md, vector = self.precipitations_presentday, type = 'node')
        if self.istemperaturescaled == 0:
            self.temperatures_reconstructed = mesh._project_3d(md, vector = self.temperatures_reconstructed, type = 'node')
        if self.isprecipscaled == 0:
            self.precipitations_reconstructed = mesh._project_3d(md, vector = self.precipitations_reconstructed, type = 'node')
        self.s0p = mesh._project_3d(md, vector = self.s0p, type = 'node')
        self.s0t = mesh._project_3d(md, vector = self.s0t, type = 'node')            

        return self
    
    def checkconsistency(self, md, solution, analyses):
        if 'MasstransportAnalysis' in analyses:
            class_utils._check_field(md, fieldname = 'smb.desfac', le = 1, scalar = True)
            class_utils._check_field(md, fieldname = 'smb.s0p', size = (md.mesh.numberofvertices, ), ge = 0, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.s0t', size = (md.mesh.numberofvertices, ), ge = 0, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.rlaps', ge = 0, scalar = True)
            class_utils._check_field(md, fieldname = 'smb.rlapslgm', ge = 0, scalar = True)
            class_utils._check_field(md, fieldname = 'smb.steps_per_step', ge = 1, scalar = True)
            class_utils._check_field(md, fieldname = 'smb.averaging', scalar = True, values = [0, 1, 2])

            if self.isd18opd:
                lent = float(np.size(self.temperatures_presentday, 1))
                lenp = float(np.size(self.precipitations_presentday, 1))
                multt = np.ceil(lent / 12.) * 12.
                multp = np.ceil(lenp / 12.) * 12.
                class_utils._check_field(md, fieldname = 'smb.temperatures_presentday', size = (md.mesh.numberofvertices, 12), allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.precipitations_presentday', size = (md.mesh.numberofvertices, 12), allow_nan = False, allow_inf = False)

                if self.istemperaturescaled == 0:
                    lent = float(np.size(self.temperatures_reconstructed, 1))
                    multt = np.ceil(lent / 12.) * 12.
                    class_utils._check_field(md, fieldname = 'smb.temperatures_reconstructed', timeseries = True, size = (md.mesh.numberofvertices + 1, multt), allow_nan = False, allow_inf = False)

                if self.isprecipscaled == 0:
                    lenp = float(np.size(self.precipitations_reconstructed, 1))
                    multp = np.ceil(lent / 12.) * 12.
                    class_utils._check_field(md, fieldname = 'smb.precipitations_reconstructed', timeseries = True, size =  (md.mesh.numberofvertices + 1, multp), allow_nan = False, allow_inf = False)

                class_utils._check_field(md, fieldname = 'smb.delta18o', singletimeseries = True, size = (2, np.nan), allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.dpermil', ge = 0, scalar = True)
                class_utils._check_field(md, fieldname = 'smb.f', ge = 0, scalar = True)

            if self.issetpddfac:
                class_utils._check_field(md, fieldname = 'smb.pddfac_snow', ge = 0, allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.pddfac_ice', ge = 0, allow_nan = False, allow_inf = False)

        class_utils._check_field(md, fieldname = 'masstransport.requested_outputs', string_list = True)

        return md
    
    # Initialise empty fields of correct dimensions
    def _initialize(self, md):
        """
        Initialise empty fields in smb.d18opdd.
        """

        if np.all(np.isnan(self.s0p)):
            self.s0p = np.zeros((md.mesh.numberofvertices))
            warnings.warn('pyissm.model.classes.smb.d18opdd: no SMBd18opdd.s0p specified -- values set as 0.')

        if np.all(np.isnan(self.s0t)):
            self.s0t = np.zeros((md.mesh.numberofvertices))
            warnings.warn('pyissm.model.classes.smb.d18opdd: no SMBd18opdd.s0t specified -- values set as 0.')

        return self
    
    # Process requested outputs, expanding 'default' to appropriate outputs
    def _process_outputs(self,
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
        execute._write_model_field(fid, prefix, name = 'md.smb.model', data = 5, format = 'Integer')

        ## Write fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'ismungsm', format = 'Boolean')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'isd18opd', format = 'Boolean')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'issetpddfac', format = 'Boolean')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'desfac', format = 'Double')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 's0p', format = 'DoubleMat', mattype = 1)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 's0t', format = 'DoubleMat', mattype = 1)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'rlaps', format = 'Double')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'rlapslgm', format = 'Double')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'Tdiff', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'sealev', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute._write_model_field(fid, prefix, name = 'md.smb.requested_outputs', data = self._process_outputs(md), format = 'StringArray')

        ## Write conditional fields
        if self.isd18opd:
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'temperatures_presentday', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'precipitations_presentday', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'istemperaturescaled', format = 'Boolean')
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'isprecipscaled', format = 'Boolean')
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'delta18o', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'dpermil', format = 'Double')
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'f', format = 'Double')

            if self.istemperaturescaled == 0:
                execute._write_model_field(fid, prefix, obj = self, fieldname = 'temperatures_reconstructed', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

            if self.isprecipscaled == 0:
                execute._write_model_field(fid, prefix, obj = self, fieldname = 'precipitations_reconstructed', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

        if self.issetpddfac:
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'pddfac_snow', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'pddfac_ice', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

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
    mesh : ISSM mesh object
        Mesh object containing mesh information. Required for initialization.
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
        Initializes the GEMB SMB parameters, requires model object for mesh information. Optionally inherits from another instance.
    __repr__(self)
        Returns a detailed string representation of the GEMB SMB parameters.
    __str__(self)
        Returns a short string identifying the class.
    _process_outputs(self, md=None, return_default_outputs=False)
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
    md.smb = pyissm.model.classes.smb.gemb(md.mesh)
    md.smb.Ta = temperature_forcing_data
    md.smb.P = precipitation_data
    md.smb.dswrf = shortwave_radiation_data
    """

    # Initialise with default parameters
    ## NOTE: md.mesh must be specified for mesh parameters.
    def __init__(self,
                 mesh = None,
                 other = None):
        self.isgraingrowth = 1
        self.isalbedo = 1
        self.isshortwave = 1
        self.isthermal = 1
        self.isaccumulation = 1
        self.ismelt = 1
        self.isdensification = 1
        self.isturbulentflux = 1
        self.isconstrainsurfaceT = 0
        self.isdeltaLWup = 0
        self.ismappedforcing = 0
        self.isprecipforcingremapped = 1
        self.iscompressedforcing = 0
        self.Ta = np.nan
        self.V = np.nan
        self.dswrf = np.nan
        self.dlwrf = np.nan
        self.P = np.nan
        self.eAir = np.nan
        self.pAir = np.nan
        self.Tmean = np.nan
        self.Vmean = 10 * np.ones((mesh.numberofelements,))
        self.C = np.nan
        self.Tz = np.nan
        self.Vz = np.nan
        self.aSnow = 0.85
        self.aValue = self.aSnow * np.ones(mesh.numberofelements,)
        self.teValue =  np.ones((mesh.numberofelements,))
        self.dulwrfValue = np.zeros((mesh.numberofelements,))
        self.mappedforcingpoint = np.nan
        self.mappedforcingelevation = np.nan
        self.lapseTaValue = -0.006
        self.lapsedlwrfValue = -0.032
        self.Dzini = 0.05 * np.ones((mesh.numberofelements, 2))
        self.Dini = 910.0 * np.ones((mesh.numberofelements, 2))
        self.Reini = 2.5 * np.ones((mesh.numberofelements, 2))
        self.Gdnini = 0.0 * np.ones((mesh.numberofelements, 2))
        self.Gspini = 0.0 * np.ones((mesh.numberofelements, 2))
        self.ECini = 0.0 * np.ones((mesh.numberofelements, ))
        self.Wini = 0.0 * np.ones((mesh.numberofelements, 2))
        self.Aini = self.aSnow * np.ones((mesh.numberofelements, 2))
        self.Adiffini = np.ones((mesh.numberofelements, 2))
        self.Tini = 273.15 * np.ones((mesh.numberofelements, 2))
        self.Sizeini = 2 * np.ones((mesh.numberofelements, ))
        self.aIdx = 1
        self.eIdx = 1
        self.tcIdx = 1
        self.swIdx = 0
        self.denIdx = 2
        self.dsnowIdx = 1
        self.zTop = 10 * np.ones((mesh.numberofelements,))
        self.dzTop = 0.05 * np.ones((mesh.numberofelements,))
        self.dzMin = self.dzTop / 2
        self.zY = 1.025 * np.ones((mesh.numberofelements,))
        self.zMax = 250 * np.ones((mesh.numberofelements,))
        self.zMin = 130  * np.ones((mesh.numberofelements,))
        self.outputFreq = 30
        self.dswdiffrf = 0.0 * np.ones(mesh.numberofelements,)
        self.szaValue = 0.0 * np.ones(mesh.numberofelements,)
        self.cotValue = 0.0 * np.ones(mesh.numberofelements,)
        self.ccsnowValue = 0.0 * np.ones(mesh.numberofelements,)
        self.cciceValue = 0.0 * np.ones(mesh.numberofelements,)
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
        s += '{}\n'.format(class_utils._field_display(self, 'isgraingrowth', 'run grain growth module (default true)'))
        s += '{}\n'.format(class_utils._field_display(self, 'isalbedo', 'run albedo module (default true)'))
        s += '{}\n'.format(class_utils._field_display(self, 'isshortwave', 'run short wave module (default true)'))
        s += '{}\n'.format(class_utils._field_display(self, 'isthermal', 'run thermal module (default true)'))
        s += '{}\n'.format(class_utils._field_display(self, 'isaccumulation', 'run accumulation module (default true)'))
        s += '{}\n'.format(class_utils._field_display(self, 'ismelt', 'run melting  module (default true)'))
        s += '{}\n'.format(class_utils._field_display(self, 'isdensification', 'run densification module (default true)'))
        s += '{}\n'.format(class_utils._field_display(self, 'isturbulentflux', 'run turbulant heat fluxes module (default true)'))
        s += '{}\n'.format(class_utils._field_display(self, 'isconstrainsurfaceT', 'constrain surface temperatures to air temperature, turn off EC and surface flux contribution to surface temperature change (default false)'))
        s += '{}\n'.format(class_utils._field_display(self, 'isdeltaLWup', 'set to true to invoke a bias in the long wave upward spatially, specified by dulwrfValue (default false)'))
        s += '{}\n'.format(class_utils._field_display(self,'ismappedforcing','set to true if forcing grid does not match model mesh, mapping specified by mappedforcingpoint (default false)'))
        s += '{}\n'.format(class_utils._field_display(self,'isprecipforcingremapped','set to true if ismappedforcing is true and precip should be downscaled from native grid (Default value is true)'))
        s += '{}\n'.format(class_utils._field_display(self,'iscompressedforcing','set to true to compress the input matrices when writing to binary (default false)'))
        s += '{}\n'.format(class_utils._field_display(self, 'Ta', '2 m air temperature, in Kelvin'))
        s += '{}\n'.format(class_utils._field_display(self, 'V', 'wind speed (m s-1)'))
        s += '{}\n'.format(class_utils._field_display(self, 'dswrf', 'downward shortwave radiation flux [W/m^2]'))
        s += '{}\n'.format(class_utils._field_display(self, 'dswdiffrf', 'downward diffusive portion of shortwave radiation flux (default to 0) [W/m^2]'))
        s += '{}\n'.format(class_utils._field_display(self, 'dlwrf', 'downward longwave radiation flux [W/m^2]'))
        s += '{}\n'.format(class_utils._field_display(self, 'P', 'precipitation [mm w.e. / m^2]'))
        s += '{}\n'.format(class_utils._field_display(self, 'eAir', 'screen level vapor pressure [Pa]'))
        s += '{}\n'.format(class_utils._field_display(self, 'pAir', 'surface pressure [Pa]'))
        s += '{}\n'.format(class_utils._field_display(self, 'Tmean', 'mean annual temperature [K]'))
        s += '{}\n'.format(class_utils._field_display(self, 'C', 'mean annual snow accumulation [kg m-2 yr-1]'))
        s += '{}\n'.format(class_utils._field_display(self, 'Vmean', 'mean annual temperature [m s-1] (default 10 m/s)'))
        s += '{}\n'.format(class_utils._field_display(self, 'Tz', 'height above ground at which temperature (T) was sampled [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'Vz', 'height above ground at which wind (V) eas sampled [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'zTop', 'depth over which grid length is constant at the top of the snopack (default 10) [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'dzTop', 'initial top vertical grid spacing (default .05) [m] '))
        s += '{}\n'.format(class_utils._field_display(self, 'dzMin', 'initial min vertical allowable grid spacing (default dzMin/2) [m] '))
        s += '{}\n'.format(class_utils._field_display(self, 'zMax', 'initial max model depth (default is min(thickness, 500)) [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'zMin', 'initial min model depth (default is min(thickness, 30)) [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'zY', 'stretch grid cells bellow top_z by a [top_dz * y ^ (cells bellow top_z)]'))
        s += '{}\n'.format(class_utils._field_display(self, 'InitDensityScaling', ['initial scaling factor multiplying the density of ice', 'which describes the density of the snowpack.']))
        s += '{}\n'.format(class_utils._field_display(self, 'ThermoDeltaTScaling', 'scaling factor to multiply the thermal diffusion timestep (delta t)'))
        s += '{}\n'.format(class_utils._field_display(self, 'outputFreq', 'output frequency in days (default is monthly, 30)'))
        s += '{}\n'.format(class_utils._field_display(self, 'adThresh', 'Apply aIdx method to all areas with densities below this value, or else apply direct input value from aValue, allowing albedo to be altered.'))
        s += '{}\n'.format(class_utils._field_display(self, 'aIdx', ['method for calculating albedo and subsurface absorption (default is 1)',
            '0: direct input from aValue parameter',
            '1: effective grain radius [Gardner & Sharp, 2009]',
            '2: effective grain radius [Brun et al., 1992; LeFebre et al., 2003], with swIdx=1, SW penetration follows grain size in 3 spectral bands (Brun et al., 1992)',
            '3: density and cloud amount [Greuell & Konzelmann, 1994]',
            '4: exponential time decay & wetness [Bougamont & Bamber, 2005]']))

        s += '{}\n'.format(class_utils._field_display(self, 'dulwrfValue', 'Specified bias to be applied to the outward long wave radiation at every element (W/m-2, +upward)'))
        s += '{}\n'.format(class_utils._field_display(self, 'teValue', 'Outward longwave radiation thermal emissivity forcing at every element (default in code is 1)'))
        s += '{}\n'.format(class_utils._field_display(self, 'teThresh', ['Apply eIdx method to all areas with effective grain radius above this value (mm),', 'or else apply direct input value from teValue, allowing emissivity to be altered.']))
        s += '{}\n'.format(class_utils._field_display(self, 'eIdx', ['method for calculating emissivity (default is 1)',
            '0: direct input from teValue parameter, no use of teThresh',
            '1: default value of 1, in areas with grain radius below teThresh',
            '2: default value of 1, in areas with grain radius below teThresh and areas of dry snow (not bare ice or wet) at the surface']))
        s += '{}\n'.format(class_utils._field_display(self, 'tcIdx', ['method for calculating thermal conductivity (default is 1)',
            '1: after Sturm et al, 1997',
            '2: after Calonne et al., 2011']))

        s += '{}\n'.format(class_utils._field_display(self,'mappedforcingpoint','Mapping of which forcing point will map to each mesh element for ismappedforcing option (integer). Size number of elements.'))
        s += '{}\n'.format(class_utils._field_display(self,'mappedforcingelevation','The elevation of each mapped forcing location (m above sea level) for ismappedforcing option. Size number of forcing points.'))
        s += '{}\n'.format(class_utils._field_display(self,'lapseTaValue','Temperature lapse rate of each mapped forcing location, if forcing has different grid and should be remapped for ismappedforcing option. (Default value is -0.006 K m-1., vector of mapping points)'))
        s += '{}\n'.format(class_utils._field_display(self,'lapsedlwrfValue','Longwave down lapse rate of each mapped forcing location, if forcing has different grid and should be remapped for ismappedforcing option. (Default value is -0.032 W m-2 m-1., vector of mapping points)'))

        # Snow properties init
        s += '{}\n'.format(class_utils._field_display(self, 'Dzini', 'Initial cell depth when restart [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'Dini', 'Initial snow density when restart [kg m-3]'))
        s += '{}\n'.format(class_utils._field_display(self, 'Reini', 'Initial grain size when restart [mm]'))
        s += '{}\n'.format(class_utils._field_display(self, 'Gdnini', 'Initial grain dricity when restart [-]'))
        s += '{}\n'.format(class_utils._field_display(self, 'Gspini', 'Initial grain sphericity when restart [-]'))
        s += '{}\n'.format(class_utils._field_display(self, 'ECini', 'Initial evaporation/condensation when restart [kg m-2]'))
        s += '{}\n'.format(class_utils._field_display(self, 'Wini', 'Initial snow water content when restart [kg m-2]'))
        s += '{}\n'.format(class_utils._field_display(self, 'Aini', 'Initial albedo when restart [-]'))
        s += '{}\n'.format(class_utils._field_display(self, 'Adiffini', 'Initial diffusive radiation albedo when restart (default to 1) [-]'))
        s += '{}\n'.format(class_utils._field_display(self, 'Tini', 'Initial snow temperature when restart [K]'))
        s += '{}\n'.format(class_utils._field_display(self, 'Sizeini', 'Initial number of layers when restart [-]'))

        # Additional albedo parameters
        s += '{}\n'.format(class_utils._field_display(self, 'aValue', 'Albedo forcing at every element'))
        if (self.aIdx in (1, 2) if isinstance(self.aIdx, int) 
            else list(self.aIdx) == [1, 2] if isinstance(self.aIdx, (list, np.ndarray)) 
            else False):
            s += '{}\n'.format(class_utils._field_display(self, 'aSnow', 'new snow albedo (0.64 - 0.89)'))
            s += '{}\n'.format(class_utils._field_display(self, 'aIce', 'albedo of ice (0.27-0.58)'))
            if self.aIdx == 1:
                s += '{}\n'.format(class_utils._field_display(self,'szaValue','Solar Zenith Angle [degree]'))
                s += '{}\n'.format(class_utils._field_display(self,'cotValue','Cloud Optical Thickness'))
                s += '{}\n'.format(class_utils._field_display(self,'ccsnowValue','concentration of light absorbing carbon for snow [ppm1]'))
                s += '{}\n'.format(class_utils._field_display(self,'cciceValue','concentration of light absorbing carbon for ice [ppm1]'))
        elif self.aIdx == 3:
            s += '{}\n'.format(class_utils._field_display(self, 'cldFrac', 'average cloud amount'))
        elif self.aIdx == 4:
            s += '{}\n'.format(class_utils._field_display(self, 't0wet', 'time scale for wet snow (15-21.9) [d]'))
            s += '{}\n'.format(class_utils._field_display(self, 't0dry', 'warm snow timescale (30) [d]'))
            s += '{}\n'.format(class_utils._field_display(self, 'K', 'time scale temperature coef. (7) [d]'))

        s += '{}\n'.format(class_utils._field_display(self, 'swIdx', 'apply all SW to top grid cell (0) or allow SW to penetrate surface (1) [default 0, if swIdx=1 and aIdx=2 function of effective radius (Brun et al., 1992) or else dependent on snow density (taken from Bassford, 2002)]'))
        s += '{}\n'.format(class_utils._field_display(self, 'denIdx', ['densification model to use (default is 2):',
            '1 = emperical model of Herron and Langway (1980)',
            '2 = semi-emperical model of Anthern et al. (2010)',
            '3 = DO NOT USE: physical model from Appix B of Anthern et al. (2010)',
            '4 = DO NOT USE: emperical model of Li and Zwally (2004)',
            '5 = DO NOT USE: modified emperical model (4) by Helsen et al. (2008)',
            '6 = Antarctica semi-emperical model of Ligtenberg et al. (2011)',
            '7 = Greenland semi-emperical model of Kuipers Munneke et al. (2015)']))
        s += '{}\n'.format(class_utils._field_display(self, 'dsnowIdx', ['model for fresh snow accumulation density (default is 1):',
            '0 = Original GEMB value, 150 kg/m^3',
            '1 = Antarctica value of fresh snow density, 350 kg/m^3',
            '2 = Greenland value of fresh snow density, 315 kg/m^3, Fausto et al. (2018)',
            '3 = Antarctica model of Kaspers et al. (2004), Make sure to set Vmean accurately',
            '4 = Greenland model of Kuipers Munneke et al. (2015)']))

        s += '{}\n'.format(class_utils._field_display(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(class_utils._field_display(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(class_utils._field_display(self, 'requested_outputs', 'additional outputs requested'))
        return s
    
    # Define class string
    def __str__(self):
        s = 'ISSM - smb.gemb Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [smb.gemb] fields to 3D
        """
        if np.shape(self.Ta)[0] == md.mesh.numberofelements or np.shape(self.Ta)[0] == md.mesh.numberofelements + 1 :
            self.Ta = mesh._project_3d(md, vector = self.Ta, type = 'element')
            self.V = mesh._project_3d(md, vector = self.V, type = 'element')
            self.dswrf = mesh._project_3d(md, vector = self.dswrf, type = 'element')
            self.dlwrf = mesh._project_3d(md, vector = self.dlwrf, type = 'element')
            self.P = mesh._project_3d(md, vector = self.P, type = 'element')
            self.eAir = mesh._project_3d(md, vector = self.eAir, type = 'element')
            self.pAir = mesh._project_3d(md, vector = self.pAir, type = 'element')

        if not np.isnan(self.Dzini):
            self.Dzini=mesh._project_3d(md,vector = self.Dzini, type = 'element')
        if not np.isnan(self.Dini):
            self.Dini=mesh._project_3d(md,vector = self.Dini, type = 'element')
        if not np.isnan(self.Reini):
            self.Reini=mesh._project_3d(md, vector = self.Reini, type = 'element')
        if not np.isnan(self.Gdnini):
            self.Gdnini=mesh._project_3d(md, vector = self.Gdnini, type = 'element')
        if not np.isnan(self.Gspini):
            self.Gspini=mesh._project_3d(md, vector = self.Gspini, type = 'element')
        if not np.isnan(self.ECini):
            self.ECini=mesh._project_3d(md, vector = self.ECini, type = 'element')
        if not np.isnan(self.Wini):
            self.Wini=mesh._project_3d(md, vector = self.Wini, type = 'element')
        if not np.isnan(self.Aini):
            self.Aini=mesh._project_3d(md, vector = self.Aini, type = 'element')
        if not np.isnan(self.Adiffini):
            self.Adiffini=mesh._project_3d(md, vector = self.Adiffini, type = 'element')
        if not np.isnan(self.Tini):
            self.Tini=mesh._project_3d(md, vector = self.Tini, type = 'element')

        if not np.isnan(self.dswdiffrf):
            self.dswdiffrf=mesh._project_3d(md, vector = self.dswdiffrf, type = 'element')
        if not np.isnan(self.szaValue):
            self.szaValue=mesh._project_3d(md, vector = self.szaValue, type = 'element')
        if not np.isnan(self.cotValue):
            self.cotValue=mesh._project_3d(md, vector = self.cotValue, type = 'element')
        if not np.isnan(self.ccsnowValue):
            self.ccsnowValue=mesh._project_3d(md, vector = self.ccsnowValue, type = 'element')
        if not np.isnan(self.cciceValue):
            self.cciceValue=mesh._project_3d(md, vector = self.cciceValue, type = 'element')

        if not np.isnan(self.aValue):
            self.aValue = mesh._project_3d(md, vector = self.aValue, type = 'element')
        if not np.isnan(self.teValue):
            self.teValue = mesh._project_3d(md, vector = self.teValue, type = 'element')
        if not np.isnan(self.mappedforcingpoint):
            self.mappedforcingpoint = mesh._project_3d(md, vector = self.mappedforcingpoint, type = 'element')
        
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):

        class_utils._check_field(md, fieldname = 'smb.isgraingrowth', values = [0, 1])
        class_utils._check_field(md, fieldname = 'smb.isalbedo', values = [0, 1])
        class_utils._check_field(md, fieldname = 'smb.isshortwave', values = [0, 1])
        class_utils._check_field(md, fieldname = 'smb.isthermal', values = [0, 1])
        class_utils._check_field(md, fieldname = 'smb.isaccumulation', values = [0, 1])
        class_utils._check_field(md, fieldname = 'smb.ismelt', values = [0, 1])
        class_utils._check_field(md, fieldname = 'smb.isdensification', values = [0, 1])
        class_utils._check_field(md, fieldname = 'smb.isturbulentflux', values = [0, 1])
        class_utils._check_field(md, fieldname = 'smb.isdeltaLWup', values = [0, 1])
        class_utils._check_field(md, fieldname = 'smb.isconstrainsurfaceT', values = [0, 1])
        class_utils._check_field(md, fieldname = 'smb.ismappedforcing', values = [0, 1])
        class_utils._check_field(md, fieldname = 'smb.isprecipforcingremapped', values = [0, 1])
        class_utils._check_field(md, fieldname = 'smb.iscompressedforcing', values = [0, 1])

        sizeta=np.shape(self.Ta)
        class_utils._check_field(md, fieldname = 'smb.Ta', mappedtimeseries = True, gt = 273-100, lt = 273+100, allow_nan = False, allow_inf = False) #-100/100 celsius min/max value
        class_utils._check_field(md, fieldname = 'smb.V', mappedtimeseries = True, ge = 0, lt = 45, size = sizeta, allow_nan = False, allow_inf = False) #max 500 km/h
        class_utils._check_field(md, fieldname = 'smb.dswrf', mappedtimeseries = True, ge = 0, le = 1400, size = sizeta, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.dswdiffrf', mappedtimeseries = True, ge = 0, le = 1400, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.dlwrf', mappedtimeseries = True, size = sizeta, ge = 0, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.P', mappedtimeseries = True, ge = 0, le = 200, size = sizeta, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.eAir', mappedtimeseries = True, size = sizeta, allow_nan = False, allow_inf = False)

        class_utils._check_field(md, fieldname = 'smb.Tmean', size = (sizeta[0]-1, ), gt = 273-100, lt = 273+100, allow_nan = False, allow_inf = False) #-100/100 celsius min/max value
        class_utils._check_field(md, fieldname = 'smb.C', size = (sizeta[0]-1, ), gt = 0, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.Vmean', size = (sizeta[0]-1, ), ge = 0, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.Tz', size =  (sizeta[0]-1, ), ge = 0, le = 5000, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.Vz', size = (sizeta[0]-1, ), ge = 0, le = 5000, allow_nan = False, allow_inf = False)

        class_utils._check_field(md, fieldname = 'smb.teValue', timeseries = True, ge = 0, le = 1, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.dulwrfValue', timeseries = True, allow_nan = False, allow_inf = False)

        if self.ismappedforcing:
            class_utils._check_field(md, fieldname = 'smb.mappedforcingpoint', size = (md.mesh.numberofelements, ), gt = 0, le = (sizeta[0]-1, ), allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.mappedforcingelevation', size = (sizeta[0]-1, ), allow_nan = False, allow_inf = False)
            if np.prod(np.shape(self.lapseTaValue)) == 1:
                warnings.warn('pyissm.model.classes.smb.gemb: smb.lapseTaValue is now a vector of mapped elements. Set to md.smb.lapseTaValue * np.ones(np.shape(md.smb.mappedforcingelevation))')
            if np.prod(np.shape(self.lapsedlwrfValue)) == 1:
                warnings.warn('pyissm.model.classes.smb.gemb: smb.lapsedlwrfValue is now a vector of mapped elements. Set to md.smb.lapsedlwrfValue * np.ones(np.shape(md.smb.mappedforcingelevation))')
            class_utils._check_field(md, fieldname = 'smb.lapseTaValue', size = (sizeta[0]-1, ), allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.lapsedlwrfValue', size = (sizeta[0]-1, ), allow_nan = False, allow_inf = False)

        class_utils._check_field(md, fieldname = 'smb.aIdx', values = [0, 1, 2, 3, 4], allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.eIdx', values = [0, 1, 2], allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.tcIdx', values = [1, 2], allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.swIdx', values = [0, 1], allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.denIdx', values = [1, 2, 3, 4, 5, 6, 7], allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.dsnowIdx', values = [0, 1, 2, 3, 4], allow_nan = False, allow_inf = False)

        class_utils._check_field(md, fieldname = 'smb.zTop', ge = 0, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.dzTop', gt = 0, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.dzMin', gt = 0, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.zY', ge = 1, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.outputFreq', gt = 0, lt = 10 * 365, allow_nan = False, allow_inf = False)  #10 years max
        class_utils._check_field(md, fieldname = 'smb.InitDensityScaling', ge = 0, le = 1, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.ThermoDeltaTScaling', ge = 0, le = 1, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.adThresh', ge = 0, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.teThresh', ge = 0, allow_nan = False, allow_inf = False)
        
        class_utils._check_field(md, fieldname = 'smb.aValue', timeseries = True, ge = 0, le = 1, allow_nan = False, allow_inf = True)
        if (self.aIdx in (1, 2) if isinstance(self.aIdx, int) 
            else list(self.aIdx) == [1, 2] if isinstance(self.aIdx, (list, np.ndarray)) 
            else False):
            class_utils._check_field(md, fieldname = 'smb.aSnow', ge = 0.64, le = 0.89, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.aIce', ge = 0.27, le = 0.58, allow_nan = False, allow_inf = False)
            if self.aIdx == 1:
                class_utils._check_field(md, fieldname = 'smb.szaValue', timeseries = True, ge = 0, le = 90, allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.cotValue', timeseries = True, ge = 0, allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.ccsnowValue', timeseries = True, ge = 0, allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.cciceValue', timeseries = True, ge = 0, allow_nan = False, allow_inf = False)
        elif self.aIdx == 3:
            class_utils._check_field(md, fieldname = 'smb.cldFrac', ge = 0, le = 1, allow_nan = False, allow_inf = False)
        elif self.aIdx == 4:
            class_utils._check_field(md, fieldname = 'smb.t0wet', ge = 15, le = 21.9, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.t0dry', ge = 30, le = 30, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.K', ge = 7, le = 7, allow_nan = False, allow_inf = False)

        # Check zTop is < local thickness
        he = np.sum(md.geometry.thickness[md.mesh.elements - 1], axis=1) / np.size(md.mesh.elements, 1)
        if np.any(he < self.zTop):
            raise IOError('SMBgemb consistency check error: zTop should be smaller than local ice thickness')
        class_utils._check_field(md, fieldname = 'smb.steps_per_step', ge = 1, scalar = True)
        class_utils._check_field(md, fieldname = 'smb.averaging', scalar = True, values = [0, 1, 2])
        class_utils._check_field(md, fieldname = 'smb.requested_outputs', string_list = True)

        return md

    # Initialise empty fields of correct dimensions
    def _initialize(self, md):
        """
        Initialise empty fields in smb.gemb.
        """

        warnings.warn('pyissm.model.classes.smb.gemb: No automatic initialisation possible for smb.gemb class. Ensure all fields are set correctly.')
        return self

    # Process requested outputs, expanding 'default' to appropriate outputs
    def _process_outputs(self,
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
        execute._write_model_field(fid, prefix, name = 'md.smb.model', data = 8, format = 'Integer')

        ## Write Boolean fields
        fieldnames = ['isgraingrowth', 'isalbedo', 'isshortwave', 'isthermal', 'isaccumulation',
                      'ismelt', 'isdensification', 'isturbulentflux', 'isconstrainsurfaceT',
                      'isdeltaLWup', 'ismappedforcing', 'isprecipforcingremapped']
        for field in fieldnames:
            execute._write_model_field(fid, prefix, obj = self, fieldname = field, format = 'Boolean')

        ## Write conditional compressed forcing fields
        if self.iscompressedforcing:
            writetype='CompressedMat'
        else:
            writetype='DoubleMat'

        fieldnames = ['Ta', 'V', 'dswrf', 'dswdiffrf', 'dlwrf', 'P', 'eAir', 'pAir']
        for field in fieldnames:
            execute._write_model_field(fid,prefix, obj = self, fieldname = field, format = writetype, mattype = 2, timeserieslength = np.shape(self.Ta)[0], yts = md.constants.yts)

        ## Write DoubleMat fields
        fieldnames = ['Tmean', 'C', 'Vmean', 'Tz', 'Vz', 'zTop', 'dzTop', 'dzMin', 'zY', 'zMax', 'zMin']
        for field in fieldnames:
            execute._write_model_field(fid, prefix, obj = self, fieldname = field, format = 'DoubleMat', mattype = 2)

        ## Write Integer fields
        fieldnames = ['aIdx', 'eIdx', 'tcIdx', 'swIdx', 'denIdx', 'dsnowIdx']
        for field in fieldnames:
            execute._write_model_field(fid, prefix, obj = self, fieldname = field, format = 'Integer')

        ## Write Double fields
        fieldnames = ['InitDensityScaling', 'ThermoDeltaTScaling', 'outputFreq', 'aSnow', 'aIce',
                      'cldFrac', 't0wet', 't0dry', 'K', 'adThresh', 'teThresh']
        for field in fieldnames:
            execute._write_model_field(fid, prefix, obj = self, fieldname = field, format = 'Double')

        ## Write DoubleMat fields
            ## mattype = 2
        fieldnames = ['aValue', 'teValue', 'dulwrfValue', 'szaValue', 'cotValue', 'ccsnowValue', 'cciceValue']
        for field in fieldnames:
            execute._write_model_field(fid, prefix, obj = self, fieldname = field, format = 'DoubleMat', mattype = 2, timeserieslength = md.mesh.numberofelements + 1, yts = md.constants.yts)

            ## mattype = 3
        fieldnames = ['Dzini', 'Dini', 'Reini', 'Gdnini', 'Gspini', 'Wini', 'Aini', 'Adiffini', 'Tini']
        for field in fieldnames:
            execute._write_model_field(fid, prefix, obj = self, fieldname = field, format = 'DoubleMat', mattype = 3)

        ## Write other fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'Sizeini', format = 'IntMat', mattype = 2)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'ECini', format = 'DoubleMat', mattype = 2)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')

        if self.ismappedforcing:
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'mappedforcingpoint', format ='IntMat', mattype = 2)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'mappedforcingelevation', format ='DoubleMat', mattype = 3)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'lapseTaValue', format ='DoubleMat', mattype = 3)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'lapsedlwrfValue', format ='DoubleMat', mattype = 3)

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

        execute._write_model_field(fid, prefix, name = 'md.smb.dt', data = dt, format = 'Double', scale = md.constants.yts)
        
        # Check if smb_dt goes evenly into transient core time step
        if (md.timestepping.time_step % dt >= 1e-10):
            raise IOError('smb_dt/dt = {}. The number of SMB time steps in one transient core time step has to be an an integer'.format(md.timestepping.time_step / dt))
        # Make sure that adaptive time step is off
        if md.timestepping.__class__.__name__ == 'timesteppingadaptive':
            raise IOError('GEMB cannot be run with adaptive timestepping.  Check class type of md.timestepping')

        execute._write_model_field(fid, prefix, name = 'md.smb.requested_outputs', data = self._process_outputs(md), format = 'StringArray')

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
    _process_outputs(self, md=None, return_default_outputs=False)
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
    md.smb = pyissm.model.classes.smb.gradients()
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
        self.requested_outputs = ['default']

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

        s += '{}\n'.format(class_utils._field_display(self, 'issmbgradients', 'is smb gradients method activated (0 or 1, default is 0)'))
        s += '{}\n'.format(class_utils._field_display(self, 'href', 'reference elevation from which deviation is used to calculate SMB adjustment in smb gradients method'))
        s += '{}\n'.format(class_utils._field_display(self, 'smbref', 'reference smb from which deviation is calculated in smb gradients method [m/yr ice equiv]'))
        s += '{}\n'.format(class_utils._field_display(self, 'b_pos', 'slope of hs - smb regression line for accumulation regime required if smb gradients is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'b_neg', 'slope of hs - smb regression line for ablation regime required if smb gradients is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(class_utils._field_display(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(class_utils._field_display(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.gradients Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [smb.gradients] fields to 3D
        """
        warnings.warn('pyissm.model.classes.smb.gradients._extrude: 3D extrusion not implemented for smb.gradients. Returning unchanged (2D) smb fields.')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        if 'MasstransportAnalysis' in analyses:
            class_utils._check_field(md, fieldname = 'smb.href', timeseries = True, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.smbref', timeseries = True, allow_nan = False, allow_inf = False)
            if np.max(np.max(np.abs(md.smb.smbref[0:-1,]))) < 1:
                print('!!! Warning: SMBgradients now expects smbref to be in m/yr ice eq. instead of mm/yr water eq.')
            class_utils._check_field(md, fieldname = 'smb.b_pos', timeseries = True, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.b_neg', timeseries = True, allow_nan = False, allow_inf = False)

        class_utils._check_field(md, fieldname = 'smb.steps_per_step', ge = 1, scalar = True)
        class_utils._check_field(md, fieldname = 'smb.averaging', scalar = True, values = [0, 1, 2])
        class_utils._check_field(md, fieldname = 'masstransport.requested_outputs', string_list = True)

        return md
    
    # Initialise empty fields of correct dimensions
    def _initialize(self, md):
        """
        Initialise empty fields in smb.gradients.
        """

        warnings.warn('pyissm.model.classes.smb.gradients: No automatic initialisation possible for smb.gradients class. Ensure all fields are set correctly.')
        return self
    
    # Process requested outputs, expanding 'default' to appropriate outputs
    def _process_outputs(self,
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
        execute._write_model_field(fid, prefix, name = 'md.smb.model', data = 6, format = 'Integer')

        execute._write_model_field(fid, prefix, obj = self, fieldname = 'href', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'smbref', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'b_pos', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'b_neg', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute._write_model_field(fid, prefix, name = 'md.smb.requested_outputs', data = self._process_outputs(md), format = 'StringArray')

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
    _process_outputs(self, md=None, return_default_outputs=False)
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
    md.smb = pyissm.model.classes.smb.gradientscomponents()
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

        s += '{}\n'.format(class_utils._field_display(self, 'issmbgradients', 'is smb gradients method activated (0 or 1, default is 0)'))
        s += '{}\n'.format(class_utils._field_display(self, 'accuref', ' reference value of the accumulation'))
        s += '{}\n'.format(class_utils._field_display(self, 'accualti', ' Altitude at which the accumulation is equal to the reference value'))
        s += '{}\n'.format(class_utils._field_display(self, 'accugrad', ' Gradient of the variation of the accumulation (0 for uniform accumulation)'))
        s += '{}\n'.format(class_utils._field_display(self, 'runoffref', ' reference value of the runoff m w.e. y-1 (temperature times ddf)'))
        s += '{}\n'.format(class_utils._field_display(self, 'runoffalti', ' Altitude at which the runoff is equal to the reference value'))
        s += '{}\n'.format(class_utils._field_display(self, 'runoffgrad', ' Gradient of the variation of the runoff (0 for uniform runoff) m w.e. m-1 y-1 (lapse rate times ddf)'))
        s += '{}\n'.format(class_utils._field_display(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(class_utils._field_display(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(class_utils._field_display(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.gradientscomponents Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [smb.gradientscomponents] fields to 3D
        """
        warnings.warn('pyissm.model.classes.smb.gradientscomponents._extrude: 3D extrusion not implemented for smb.gradientscomponents. Returning unchanged (2D) smb fields.')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        if 'MasstransportAnalysis' in analyses:
            class_utils._check_field(md, fieldname = 'smb.accualti', scalar = True, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.accuref', singletimeseries = True, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.accugrad', singletimeseries = True, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.runoffalti', scalar = True, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.runoffref', singletimeseries = True, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.runoffgrad', singletimeseries = True, allow_nan = False, allow_inf = False)

        class_utils._check_field(md, fieldname = 'smb.steps_per_step', ge = 1, scalar = True)
        class_utils._check_field(md, fieldname = 'smb.averaging', scalar = True, values = [0, 1, 2])
        class_utils._check_field(md, fieldname = 'masstransport.requested_outputs', string_list = True)
        
        return md
    
    # Initialise empty fields of correct dimensions
    def _initialize(self, md):
        """
        Initialise empty fields in smb.gradientscomponents.
        """

        warnings.warn('pyissm.model.classes.smb.gradientscomponents: No automatic initialisation possible for smb.gradientscomponents class. Ensure all fields are set correctly.')
        return self
    
    # Process requested outputs, expanding 'default' to appropriate outputs
    def _process_outputs(self,
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
        execute._write_model_field(fid, prefix, name = 'md.smb.model', data = 11, format = 'Integer')

        ## Write DoubleMat fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'accuref', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts, scale = 1. / md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'accugrad', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts, scale = 1. / md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'runoffref', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts, scale = 1. / md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'runoffgrad', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts, scale = 1. / md.constants.yts)
        
        ## Write other fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'accualti', format = 'Double')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'runoffalti', format = 'Double')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute._write_model_field(fid, prefix, name = 'md.smb.requested_outputs', data = self._process_outputs(md), format = 'StringArray')

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
    _process_outputs(self, md=None, return_default_outputs=False)
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
    md.smb = pyissm.model.classes.smb.gradientsela()
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
        s += '{}\n'.format(class_utils._field_display(self, 'ela', ' equilibrium line altitude from which deviation is used to calculate smb using the smb gradients ela method [m a.s.l.]'))
        s += '{}\n'.format(class_utils._field_display(self, 'b_pos', ' vertical smb gradient (dB/dz) above ela'))
        s += '{}\n'.format(class_utils._field_display(self, 'b_neg', ' vertical smb gradient (dB/dz) below ela'))
        s += '{}\n'.format(class_utils._field_display(self, 'b_max', ' upper cap on smb rate, default: 9999 (no cap) [m ice eq./yr]'))
        s += '{}\n'.format(class_utils._field_display(self, 'b_min', ' lower cap on smb rate, default: -9999 (no cap) [m ice eq./yr]'))
        s += '{}\n'.format(class_utils._field_display(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(class_utils._field_display(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(class_utils._field_display(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.gradientsela Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [smb.gradientsela] fields to 3D
        """
        warnings.warn('pyissm.model.classes.smb.gradientsela._extrude: 3D extrusion not implemented for smb.gradientsela. Returning unchanged (2D) smb fields.')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        if 'MasstransportAnalysis' in analyses:
            class_utils._check_field(md, fieldname = 'smb.ela', timeseries = True, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.b_pos', timeseries = True, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.b_neg', timeseries = True, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.b_max', timeseries = True, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.b_min', timeseries = True, allow_nan = False, allow_inf = False)

        class_utils._check_field(md, fieldname = 'smb.steps_per_step', ge = 1, scalar = True)
        class_utils._check_field(md, fieldname = 'smb.averaging', scalar = True, values = [0, 1, 2])
        class_utils._check_field(md, fieldname = 'smb.requested_outputs', string_list = True)

        return md
    
    # Initialise empty fields of correct dimensions
    def _initialize(self, md):   
        """
        Initialise empty fields in smb.gradientsela.
        """

        warnings.warn('pyissm.model.classes.smb.gradientsela: No automatic initialisation possible for smb.gradientsela class. Ensure all fields are set correctly.')
        return self
    
    # Process requested outputs, expanding 'default' to appropriate outputs
    def _process_outputs(self,
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
        execute._write_model_field(fid, prefix, name = 'md.smb.model', data = 9, format = 'Integer')

        ## Write DoubleMat fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'ela', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'b_pos', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'b_neg', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'b_max', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'b_min', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        
        ## Write other fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute._write_model_field(fid, prefix, name = 'md.smb.requested_outputs', data = self._process_outputs(md), format = 'StringArray')

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
    _process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.smb = pyissm.model.classes.smb.henning()
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

        s += '{}\n'.format(class_utils._field_display(self, 'smbref', 'reference smb from which deviation is calculated [m/yr ice eq]'))
        s += '{}\n'.format(class_utils._field_display(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(class_utils._field_display(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(class_utils._field_display(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.henning Class'
        return s
        
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [smb.henning] fields to 3D
        """
        self.smbref = mesh._project_3d(md, vector = self.smbref, type = 'node')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        if solution == 'TransientSolution' and not md.transient.issmb:
            return
        if 'MasstransportAnalysis' in analyses:
            class_utils._check_field(md, fieldname = 'smb.mass_balance', timeseries = True, allow_nan = False, allow_inf = False)
        if 'BalancethicknessAnalysis' in analyses:
            class_utils._check_field(md, fieldname = 'smb.mass_balance', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'smb.steps_per_step', ge = 1, scalar = True)
        class_utils._check_field(md, fieldname = 'smb.averaging', scalar = True, values = [0, 1, 2])
        class_utils._check_field(md, fieldname = 'smb.requested_outputs', string_list = True)

        return md
    
    # Initialise empty fields of correct dimensions
    def _initialize(self, md):
        """
        Initialise empty fields in smb.henning.
        """

        if np.all(np.isnan(self.smbref)):
            self.smbref = np.zeros((md.mesh.numberofvertices))
            warnings.warn('pyissm.model.classes.smb.henning: no smb.smbref specified -- values set as 0.')

        return self

    # Process requested outputs, expanding 'default' to appropriate outputs
    def _process_outputs(self,
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
        execute._write_model_field(fid, prefix, name = 'md.smb.model', data = 7, format = 'Integer')

        ## Write fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'smbref', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute._write_model_field(fid, prefix, name = 'md.smb.requested_outputs', data = self._process_outputs(md), format = 'StringArray')

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
    _process_outputs(self, md=None, return_default_outputs=False)
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
    md.smb = pyissm.model.classes.smb.meltcomponents()
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

        s += '{}\n'.format(class_utils._field_display(self, 'accumulation', 'accumulated snow [m/yr ice eq]'))
        s += '{}\n'.format(class_utils._field_display(self, 'evaporation', 'mount of ice lost to evaporative processes [m/yr ice eq]'))
        s += '{}\n'.format(class_utils._field_display(self, 'melt', 'amount of ice melt in the ice column [m/yr ice eq]'))
        s += '{}\n'.format(class_utils._field_display(self, 'refreeze', 'amount of ice melt refrozen in the ice column [m/yr ice eq]'))
        s += '{}\n'.format(class_utils._field_display(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(class_utils._field_display(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(class_utils._field_display(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.meltcomponents Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [smb.meltcomponents] fields to 3D
        """
        self.accumulation = mesh._project_3d(md, vector = self.accumulation, type = 'node')
        self.evaporation = mesh._project_3d(md, vector = self.evaporation, type = 'node')
        self.melt = mesh._project_3d(md, vector = self.melt, type = 'node')
        self.refreeze = mesh._project_3d(md, vector = self.refreeze, type = 'node')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        if 'MasstransportAnalysis' in analyses:
            class_utils._check_field(md, fieldname = 'smb.accumulation', timeseries = True, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.evaporation', timeseries = True, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.refreeze', timeseries = True, allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.melt', timeseries = True, allow_nan = False, allow_inf = False)

        if 'BalancethicknessAnalysis' in analyses:
            class_utils._check_field(md, fieldname = 'smb.accumulation', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.evaporation', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.refreeze', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.melt', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)

        class_utils._check_field(md, fieldname = 'smb.steps_per_step', ge = 1, scalar = True)
        class_utils._check_field(md, fieldname = 'smb.averaging', scalar = True, values = [0, 1, 2])
        class_utils._check_field(md, fieldname = 'smb.requested_outputs', string_list = 1)

        return md
    
    # Initialise empty fields of correct dimensions
    def _initialize(self, md):
        """
        Initialise empty fields in smb.meltcomponents.
        """
        
        if np.all(np.isnan(self.accumulation)):
            self.accumulation = np.zeros((md.mesh.numberofvertices))
            warnings.warn('pyissm.model.classes.smb.meltcomponents: no smb.accumulation specified -- values set as 0.')

        if np.all(np.isnan(self.evaporation)):
            self.evaporation = np.zeros((md.mesh.numberofvertices))
            warnings.warn('pyissm.model.classes.smb.meltcomponents: no smb.evaporation specified -- values set as 0.')

        if np.all(np.isnan(self.refreeze)):
            self.refreeze = np.zeros((md.mesh.numberofvertices))
            warnings.warn('pyissm.model.classes.smb.meltcomponents: no smb.refreeze specified -- values set as 0.')

        if np.all(np.isnan(self.melt)):
            self.melt = np.zeros((md.mesh.numberofvertices))
            warnings.warn('pyissm.model.classes.smb.meltcomponents: no smb.melt specified -- values set as 0.')

        return self

    # Process requested outputs, expanding 'default' to appropriate outputs
    def _process_outputs(self,
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
        execute._write_model_field(fid, prefix, name = 'md.smb.model', data = 3, format = 'Integer')

        ## Write fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'accumulation', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'evaporation', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'melt', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'refreeze', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute._write_model_field(fid, prefix, name = 'md.smb.requested_outputs', data = self._process_outputs(md), format = 'StringArray')

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
    _process_outputs(self, md=None, return_default_outputs=False)
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
    md.smb = pyissm.model.classes.smb.pdd()
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
        self.requested_outputs = ['default']

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'

        s += '{}\n'.format(class_utils._field_display(self, 'isdelta18o', 'is temperature and precipitation delta18o parametrisation activated (0 or 1, default is 0)'))
        s += '{}\n'.format(class_utils._field_display(self, 'ismungsm', 'is temperature and precipitation mungsm parametrisation activated (0 or 1, default is 0)'))
        s += '{}\n'.format(class_utils._field_display(self, 'issetpddfac', 'is user passing in defined pdd factors (0 or 1, default is 0)'))
        s += '{}\n'.format(class_utils._field_display(self, 'desfac', 'desertification elevation factor (between 0 and 1, default is 0.5) [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 's0p', 'should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 's0t', 'should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'rlaps', 'present day lapse rate [degree/km]'))
        s += '{}\n'.format(class_utils._field_display(self, 'rlapslgm', 'LGM lapse rate [degree/km]'))
        s += '{}\n'.format(class_utils._field_display(self, 'monthlytemperatures', 'monthly surface temperatures [K], required if pdd is activated and delta18o not activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'precipitation', 'monthly surface precipitation [m/yr water eq], required if pdd is activated and delta18o or mungsm not activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'delta18o', 'delta18o [per mil], required if pdd is activated and delta18o activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'delta18o_surface', 'surface elevation of the delta18o site, required if pdd is activated and delta18o activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'temperatures_presentday', 'monthly present day surface temperatures [K], required if delta18o/mungsm is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'temperatures_lgm', 'monthly LGM surface temperatures [K], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'precipitations_presentday', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'precipitations_lgm', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'Tdiff', 'time interpolation parameter for temperature, 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'sealev', 'sea level [m], 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'temperatures_presentday', 'monthly present day surface temperatures [K], required if delta18o/mungsm is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'temperatures_lgm', 'monthly LGM surface temperatures [K], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'precipitations_presentday', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'precipitations_lgm', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'Pfac', 'time interpolation parameter for precipitation, 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'Tdiff', 'time interpolation parameter for temperature, 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'sealev', 'sea level [m], 1D(year), required if mungsm is activated'))
        s += '{}\n'.format(class_utils._field_display(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(class_utils._field_display(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(class_utils._field_display(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.pdd Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [smb.pdd] fields to 3D
        """
        if not (self.isdelta18o and self.ismungsm):
            self.precipitation = mesh._project_3d(md, vector = self.precipitation, type = 'node')
            self.monthlytemperatures = mesh._project_3d(md, vector = self.monthlytemperatures, type = 'node')

        if self.isdelta18o:
            self.temperatures_lgm = mesh._project_3d(md, vector = self.temperatures_lgm, type = 'node')
            self.temperatures_presentday = mesh._project_3d(md, vector = self.temperatures_presentday, type = 'node')
            self.precipitations_presentday = mesh._project_3d(md, vector = self.precipitations_presentday, type = 'node')
            self.precipitations_lgm = mesh._project_3d(md, vector = self.precipitations_lgm, type = 'node')

        if self.ismungsm:
            self.temperatures_lgm = mesh._project_3d(md, vector = self.temperatures_lgm, type = 'node')
            self.temperatures_presentday = mesh._project_3d(md, vector = self.temperatures_presentday, type = 'node')
            self.precipitations_presentday = mesh._project_3d(md, vector = self.precipitations_presentday, type = 'node')
            self.precipitations_lgm = mesh._project_3d(md, vector = self.precipitations_lgm, type = 'node')

        if self.issetpddfac:
            self.pddfac_snow = mesh._project_3d(md, vector = self.pddfac_snow, type = 'node')
        if self.issetpddfac:
            self.pddfac_ice = mesh._project_3d(md, vector = self.pddfac_ice, type = 'node')
        self.s0p = mesh._project_3d(md, vector = self.s0p, type = 'node')
        self.s0t = mesh._project_3d(md, vector = self.s0t, type = 'node')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):

        if 'MasstransportAnalysis' in analyses:
            class_utils._check_field(md, fieldname = 'smb.desfac', le = 1, scalar = True)
            class_utils._check_field(md, fieldname = 'smb.s0p', ge = 0, size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.s0t', ge = 0, size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.rlaps', ge = 0, scalar = True)
            class_utils._check_field(md, fieldname = 'smb.rlapslgm', ge = 0, scalar = True)

            if (self.isdelta18o == 0 and self.ismungsm == 0):
                class_utils._check_field(md, fieldname = 'smb.monthlytemperatures', size = (md.mesh.numberofvertices, 12), allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.precipitation', size = (md.mesh.numberofvertices, 12), allow_nan = False, allow_inf = False)
            elif self.isdelta18o:
                class_utils._check_field(md, fieldname = 'smb.delta18o', size = (2, np.nan, ), singletimeseries = True, allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.delta18o_surface', size = (2, np.nan, ), singletimeseries = True, allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.temperatures_presentday', size = (md.mesh.numberofvertices, 12), allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.temperatures_lgm', size = (md.mesh.numberofvertices, 12), allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.precipitations_presentday', size = (md.mesh.numberofvertices, 12), allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.precipitations_lgm', size = (md.mesh.numberofvertices, 12), allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.Tdiff', size = (2, np.nan), singletimeseries = True, allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.sealev', size = (2, np.nan), singletimeseries = True, allow_nan = False, allow_inf = False)
            elif self.ismungsm:
                class_utils._check_field(md, fieldname = 'smb.temperatures_presentday', size = (md.mesh.numberofvertices, 12), allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.temperatures_lgm', size = (md.mesh.numberofvertices, 12), allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.precipitations_presentday', size = (md.mesh.numberofvertices, 12), allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.precipitations_lgm', size = (md.mesh.numberofvertices, 12), allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.Pfac', size = (2, np.nan), singletimeseries = True, allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.Tdiff', size = (2, np.nan), singletimeseries = True, allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.sealev', size = (2, np.nan), singletimeseries = True, allow_nan = False, allow_inf = False)

            if self.issetpddfac:
                class_utils._check_field(md, fieldname = 'smb.pddfac_snow', ge = 0, allow_nan = False, allow_inf = False)
                class_utils._check_field(md, fieldname = 'smb.pddfac_ice', ge = 0, allow_nan = False, allow_inf = False)

        class_utils._check_field(md, fieldname = 'smb.steps_per_step', ge = 1, scalar = True)
        class_utils._check_field(md, fieldname = 'smb.averaging', scalar = True, values = [0, 1, 2])
        class_utils._check_field(md, fieldname = 'masstransport.requested_outputs', string_list = 1)

        return md
    
    # Initialise empty fields of correct dimensions
    def _initialize(self, md):
        """
        Initialise empty fields in smb.pdd.
        """
        
        if np.all(np.isnan(self.s0p)):
            self.s0p = np.zeros((md.mesh.numberofvertices))
            warnings.warn('pyissm.model.classes.smb.pdd: no SMBpdd.s0p specified -- values set as 0.')

        if np.all(np.isnan(self.s0t)):
            self.s0t = np.zeros((md.mesh.numberofvertices))
            warnings.warn('pyissm.model.classes.smb.pdd: no SMBpdd.s0t specified -- values set as 0.')

        return self

    # Process requested outputs, expanding 'default' to appropriate outputs
    def _process_outputs(self,
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
        execute._write_model_field(fid, prefix, name = 'md.smb.model', data = 4, format = 'Integer')

        ## Write fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'isdelta18o', format = 'Boolean')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'ismungsm', format = 'Boolean')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'issetpddfac', format = 'Boolean')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'desfac', format = 'Double')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 's0p', format = 'DoubleMat', mattype = 1)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 's0t', format = 'DoubleMat', mattype = 1)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'rlaps', format = 'Double')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'rlapslgm', format = 'Double')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute._write_model_field(fid, prefix, name = 'md.smb.requested_outputs', data = self._process_outputs(md), format = 'StringArray')

        ## Write conditional fields
        if (self.isdelta18o == 0 and self.ismungsm == 0):
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'monthlytemperatures', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'precipitation', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        elif self.isdelta18o:
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'temperatures_presentday', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'temperatures_lgm', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'precipitations_presentday', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'precipitations_lgm', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'delta18o_surface', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'delta18o', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'Tdiff', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'sealev', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)
        elif self.ismungsm:
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'temperatures_presentday', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'temperatures_lgm', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'precipitations_presentday', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'precipitations_lgm', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'Pfac', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'Tdiff', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'sealev', format = 'DoubleMat', mattype = 1, timeserieslength = 2, yts = md.constants.yts)

        if self.issetpddfac:
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'pddfac_snow', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'pddfac_ice', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

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
    pdd_fac_ice : float, default=7.28
         Pdd factor for ice for all the domain (mm ice equiv/day/decgree C)
    pdd_fac_snow : float, default=2.73
         Pdd factor for snow for all the domain (mm ice equiv/day/decgree C)
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
    _process_outputs(self, md=None, return_default_outputs=False)
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
    md.smb = pyissm.model.classes.smb.pddSicopolis()
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
        self.pdd_fac_ice = 7.28
        self.pdd_fac_snow = 2.73
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = ['default']

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'
        s += '   SICOPOLIS PDD scheme (Calov & Greve, 2005):\n'
        s += '{}\n'.format(class_utils._field_display(self, 'monthlytemperatures', 'monthly surface temperatures [K]'))
        s += '{}\n'.format(class_utils._field_display(self, 'precipitation', 'monthly surface precipitation [m/yr water eq]'))
        s += '{}\n'.format(class_utils._field_display(self, 'temperature_anomaly', 'anomaly to monthly reference temperature (additive [K])'))
        s += '{}\n'.format(class_utils._field_display(self, 'precipitation_anomaly', 'anomaly to monthly precipitation (multiplicative, e.g. q = q0*exp(0.070458*DeltaT) after Huybrechts (2002)) [no unit])'))
        s += '{}\n'.format(class_utils._field_display(self, 'smb_corr', 'correction of smb after PDD call [m/a]'))
        s += '{}\n'.format(class_utils._field_display(self, 's0p', 'should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 's0t', 'should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'rlaps', 'present day lapse rate (default is 7.4 degree/km)'))
        s += '{}\n'.format(class_utils._field_display(self, 'desfac', 'desertification elevation factor (default is -log(2.0)/1000)'))
        s += '{}\n'.format(class_utils._field_display(self, 'isfirnwarming', 'is firnwarming (Reeh 1991) activated (0 or 1, default is 1)'))
        s += '{}\n'.format(class_utils._field_display(self, 'pdd_fac_ice', 'Pdd factor for ice for all the domain [mm ice equiv/day/decgree C]'))
        s += '{}\n'.format(class_utils._field_display(self, 'pdd_fac_snow', 'Pdd factor for snow for all the domain [mm ice equiv/day/decgree C]'))
        s += '{}\n'.format(class_utils._field_display(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(class_utils._field_display(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(class_utils._field_display(self, 'requested_outputs', 'additional outputs requested (TemperaturePDD, SmbAccumulation, SmbMelt)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.pddSicopolis Class'
        return s

    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [smb.pddSicopolis] fields to 3D
        """
        self.precipitation = mesh._project_3d(md, vector = self.precipitation, type = 'node')
        self.monthlytemperatures = mesh._project_3d(md, vector = self.monthlytemperatures, type = 'node')
        self.temperature_anomaly = mesh._project_3d(md, vector = self.temperature_anomaly, type = 'node')
        self.precipitation_anomaly = mesh._project_3d(md, vector = self.precipitation_anomaly, type = 'node')
        self.smb_corr = mesh._project_3d(md, vector = self.smb_corr, type = 'node')
        self.s0p = mesh._project_3d(md, vector = self.s0p, type = 'node')
        self.s0t = mesh._project_3d(md, vector = self.s0t, type = 'node')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        if solution == 'TransientSolution' and not md.transient.issmb:
            return
        if 'MasstransportAnalysis' in analyses:
            class_utils._check_field(md, fieldname = 'smb.desfac', le = 1, scalar = True)
            class_utils._check_field(md, fieldname = 'smb.s0p', ge = 0, size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.s0t', ge = 0, size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.rlaps', ge = 0, scalar = True)
            class_utils._check_field(md, fieldname = 'smb.monthlytemperatures', size = (md.mesh.numberofvertices, 12), allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.precipitation', size = (md.mesh.numberofvertices, 12), allow_nan = False, allow_inf = False)
            class_utils._check_field(md, fieldname = 'smb.pdd_fac_ice', gt = 0, scalar = True)
            class_utils._check_field(md, fieldname = 'smb.pdd_fac_snow', gt = 0, scalar = True)
        class_utils._check_field(md, fieldname = 'smb.steps_per_step', ge = 1, scalar = True)
        class_utils._check_field(md, fieldname = 'smb.averaging', scalar = True, values = [0, 1, 2])
        class_utils._check_field(md, fieldname = 'smb.requested_outputs', string_list = 1)
        
        return md
    
    # Initialise empty fields of correct dimensions
    def _initialize(self, md):
        """
        Initialise empty fields in smb.pddSicopolis.
        """

        if np.isnan(self.s0p):
            self.s0p = np.zeros((md.mesh.numberofvertices, ))
            warnings.warn('pyissm.model.classes.smb.pddSicopolis: no SMBpddSicopolis.s0p specified -- values set as 0.')

        if np.isnan(self.s0t):
            self.s0t = np.zeros((md.mesh.numberofvertices, ))
            warnings.warn('pyissm.model.classes.smb.pddSicopolis: no SMBpddSicopolis.s0t specified -- values set as 0.')

        if np.isnan(self.temperature_anomaly):
            self.temperature_anomaly = np.zeros((md.mesh.numberofvertices, ))
            warnings.warn('pyissm.model.classes.smb.pddSicopolis: no SMBpddSicopolis.temperature_anomaly specified -- values set as 0.')

        if np.isnan(self.precipitation_anomaly):
            self.precipitation_anomaly = np.ones((md.mesh.numberofvertices, ))
            warnings.warn('pyissm.model.classes.smb.pddSicopolis: no SMBpddSicopolis.precipitation_anomaly specified -- values set as 1.')

        if np.isnan(self.smb_corr):
            self.smb_corr = np.zeros((md.mesh.numberofvertices, ))
            warnings.warn('pyissm.model.classes.smb.pddSicopolis: no SMBpddSicopolis.smb_corr specified -- values set as 0.')

        return self

    # Process requested outputs, expanding 'default' to appropriate outputs
    def _process_outputs(self,
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
        execute._write_model_field(fid, prefix, name = 'md.smb.model', data = 10, format = 'Integer')

        ## Write fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'isfirnwarming', format = 'Boolean')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'desfac', format = 'Double')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 's0p', format = 'DoubleMat', mattype = 1)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 's0t', format = 'DoubleMat', mattype = 1)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'rlaps', format = 'Double')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'pdd_fac_ice', format = 'Double')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'pdd_fac_snow', format = 'Double')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'rlaps', format = 'Double')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'monthlytemperatures', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'precipitation', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'temperature_anomaly', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'precipitation_anomaly', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'smb_corr', format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'steps_per_step', format = 'Integer')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'averaging', format = 'Integer')
        execute._write_model_field(fid, prefix, name = 'md.smb.requested_outputs', data = self._process_outputs(md), format = 'StringArray')

## ------------------------------------------------------
## smb.semic
## ------------------------------------------------------
@class_registry.register_class
class semic(class_registry.manage_state):
    
    # Initialise with default parameters
    def __init__(self, other = None):
        self.dailysnowfall = np.nan
        self.dailyrainfall = np.nan
        self.dailydsradiation = np.nan
        self.dailydlradiation = np.nan
        self.dailywindspeed = np.nan
        self.dailypressure = np.nan
        self.dailyairdensity = np.nan
        self.dailyairhumidity = np.nan
        self.dailytemperature = np.nan
        self.Tamp = np.nan
        self.mask = np.nan
        self.hice = np.nan
        self.hsnow = np.nan
        self.desfac = -np.log(2.0) / 1000
        self.desfacElevation = 2000
        self.rlaps = 7.4
        self.rdl = 29
        self.s0gcm = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.hcrit = 0.028
        self.rcrit = 0.85
        self.albedo = 0
        self.albedo_snow = 0
        self.albedo_scheme = 0
        self.alb_smax = 0.79
        self.alb_smin = 0.6
        self.albi = 0.41
        self.albl = 0.07
        self.tmin = 263.15
        self.tmax = 273.15
        self.mcrit = 6e-8
        self.tau_a = 0.008
        self.tau_f = 0.24
        self.wcrit = 15.0
        self.tmid = 273.35
        self.afac = -0.18
        self.ismethod = 0
        self.isdesertification = 1
        self.isLWDcorrect = 1
        self.requested_outputs = ['default']

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surface forcings parameters:\n'
        s += '   Interface for coupling GCM data to the energy balance model SEMIC (Krapp et al (2017) https://doi.org/10.5194/tc-11-1519-2017).\n'
        s += '   The implemented coupling uses daily mean GCM input to calculate yearly mean smb, accumulation, ablation, and surface temperature.\n'
        s += '   smb and temperatures are updated every year\n'
        s += '\n   SEMIC parameters:\n'
        s += '{}\n'.format((self, 'dailysnowfall', 'daily surface dailysnowfall [m/s]'))
        s += '{}\n'.format(class_utils._field_display(self, 'dailyrainfall', 'daily surface dailyrainfall [m/s]'))
        s += '{}\n'.format(class_utils._field_display(self, 'dailydsradiation', 'daily downwelling shortwave radiation [W/m2]'))
        s += '{}\n'.format(class_utils._field_display(self, 'dailydlradiation', 'daily downwelling longwave radiation [W/m2]'))
        s += '{}\n'.format(class_utils._field_display(self, 'dailywindspeed', 'daily surface wind speed [m/s]'))
        s += '{}\n'.format(class_utils._field_display(self, 'dailypressure', 'daily surface pressure [Pa]'))
        s += '{}\n'.format(class_utils._field_display(self, 'dailyairdensity', 'daily air density [kg/m3]'))
        s += '{}\n'.format(class_utils._field_display(self, 'dailyairhumidity', 'daily air specific humidity [kg/kg]'))
        s += '{}\n'.format(class_utils._field_display(self, 'rlaps', 'present day lapse rate (default is 7.4 [degree/km]; Erokhina et al. 2017)'))
        s += '{}\n'.format(class_utils._field_display(self, 'desfac', 'desertification elevation factor (default is -log(2.0)/1000 [1/m]; Vizcaino et al. 2010)'))
        s += '{}\n'.format(class_utils._field_display(self, 'rdl', 'longwave downward radiation decrease (default is 29 [W/m^2/km]; Marty et al. 2002)'))
        s += '{}\n'.format(class_utils._field_display(self, 's0gcm', 'GCM reference elevation; (default is 0) [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'ismethod','method for calculating SMB with SEMIC. Default version of SEMIC is really slow. 0: steady, 1: transient (default: 0)'))
        if self.ismethod: # transient mode
            s += '{}\n'.format(class_utils._field_display(self,'desfacElevation','desertification elevation (default is 2000 m; Vizcaino et al. 2010)'))
            s += '{}\n'.format(class_utils._field_display(self,'Tamp','amplitude of diurnal cycle [K]'))
            s += '{}\n'.format(class_utils._field_display(self,'albedo','initial albedo [no unit]'))
            s += '{}\n'.format(class_utils._field_display(self,'albedo_snow','initial albedo for snow [no unit]'))
            s += '{}\n'.format(class_utils._field_display(self,'hice','initial thickness of ice [unit: m]'))
            s += '{}\n'.format(class_utils._field_display(self,'hsnow','initial thickness of snow [unit: m]'))
            s += '{}\n'.format(class_utils._field_display(self,'mask','masking for albedo. 0: ocean, 1: land, 2: ice (default: 2)'))
            s += '{}\n'.format(class_utils._field_display(self,'qmr','initial net energy difference between melt and refreeze in SEMIC [unit: W m^{-2}]. This variable can be set with zeros because net energy difference between melt and refreeze is dissipated fast.'))
            s += '{}\n'.format(class_utils._field_display(self,'hcrit','critical snow height for albedo [unit: m]'))
            s += '{}\n'.format(class_utils._field_display(self,'rcrit','critical refreezing height for albedo [no unit]'))

            s += '\nSEMIC albedo parameters.\n'
            s += '{}\n'.format(class_utils._field_display(self,'albedo_scheme','albedo scheme for SEMIC. 0: none, 1: slater, 2: denby, 3: isba, 4: alex (default is 0)'))
            s += '{}\n'.format(class_utils._field_display(self,'alb_smax','maximum snow albedo (default: 0.79)'))
            s += '{}\n'.format(class_utils._field_display(self,'alb_smin','minimum snow albedo (default: 0.6)'))
            s += '{}\n'.format(class_utils._field_display(self,'albi','background albedo for bare ice (default: 0.41)'))
            s += '{}\n'.format(class_utils._field_display(self,'albl','background albedo for bare land (default: 0.07)'))
            
            s += '{}\n'.format(class_utils._field_display(self,'isdesertification','enable or disable desertification of Vizcaino et al. (2010). 0: off, 1: on (default: 1)'))
            s += '{}\n'.format(class_utils._field_display(self,'isLWDcorrect','enable or disable downward longwave correction of Marty et al. (2002). 0: off, 1: on (default: 1)'))
        # albedo_scheme - 0: none, 1: slater, 2: isba, 3: denby, 4: alex.
        if self.albedo_scheme == 0:
            s += '\n\tSEMIC snow albedo parameter of None.\n'
            s += '\t   albedo of snow is updated from albedo snow max (alb_smax).\n'
            s += '\t   alb_snow = abl_smax \n '
        elif self.albedo_scheme == 1:
            s += '\n\tSEMIC snow albedo parameters of Slater et al, (1998).\n'
            s += '\t   alb = alb_smax - (alb_smax - alb_smin)*tm^(3.0)\n'
            s += '\t   tm  = 1 (tsurf > 273.15 K)\n'
            s += '\t         tm = f*(tsurf-tmin) (tmin <= tsurf < 273.15)\n'
            s += '\t         0 (tsurf < tmin)\n'
            s += '\t   f = 1/(273.15-tmin)\n'
            s += '{}\n'.format(class_utils._field_display(self, 'tmin', 'minimum temperature for which albedo decline become effective. (default: 263.15 K)[unit: K])'))
            s += '{}\n'.format(class_utils._field_display(self, 'tmax', 'maxmium temperature for which albedo decline become effective. This value should be fixed. (default: 273.15 K)[unit: K])'))
        elif self.albedo_scheme == 2:
            s += '\n\tSEMIC snow albedo parameters of Denby et al. (2002 Tellus)\n'
            s += '{}\n'.format(class_utils._field_display(self,'mcrit','critical melt rate (defaut: 6e-8) [unit: m/sec]'))
        elif self.albedo_scheme == 3:
            s += '\n\tSEMIC snow albedo parameters of ISB (Douville et al., 1995).\n'
            s += '{}\n'.format(class_utils._field_display(self, 'mcrit', 'critical melt rate (default: 6e-8) [unit: m/sec]'))
            s += '{}\n'.format(class_utils._field_display(self, 'wcrit', 'critical liquid water content (default: 15) [unit: kg/m2]'))
            s += '{}\n'.format(class_utils._field_display(self, 'tau_a', 'dry albedo decline [unit: 1/day]'))
            s += '{}\n'.format(class_utils._field_display(self, 'tau_f', 'wet albedo decline [unit: 1/day]'))
            s += '\n\tReference'
            s += '\tDouville, H., Royer, J.-F., and Mahfouf, J.-F.: A new snow parameterization for the Météo-France climate model. Part I: validation in stand-alone experiments, Climate Dynamics, 12, 21–35, https://doi.org/10.1007/s003820050092, 1995.'
        elif self.albedo_scheme == 4:
            s += '\n\tSEMIC snow albedo parameters of Alex.?\n'
            s += '{}\n'.format(class_utils._field_display(self,'afac','[unit: ?]'))
            s += '{}\n'.format(class_utils._field_display(self,'tmid','[unit: ?]'))
        else:
            raise Exception('ERROR: {} is not supported albedo scheme.'.format(self.albedo_scheme))

        s += '{}\n'.format(class_utils._field_display(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(class_utils._field_display(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(class_utils._field_display(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - smb.semic Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [smb.semic] fields to 3D
        """
        self.dailysnowfall = mesh.project3d(md, vector = self.dailysnowfall, type = 'node')
        self.dailyrainfall = mesh.project3d(md, vector = self.dailyrainfall, type = 'node')
        self.dailydsradiation = mesh.project3d(md, vector = self.dailydsradiation, type = 'node')
        self.dailydlradiation = mesh.project3d(md, vector = self.dailydlradiation, type = 'node')
        self.dailywindspeed = mesh.project3d(md, vector = self.dailywindspeed, type = 'node')
        self.dailypressure = mesh.project3d(md, vector = self.dailypressure, type = 'node')
        self.dailyairdensity = mesh.project3d(md, vector = self.dailyairdensity, type = 'node')
        self.dailyairhumidity = mesh.project3d(md, vector = self.dailyairhumidity, type = 'node')
        self.dailytemperature = mesh.project3d(md, vector = self.dailytemperature, type = 'node')
        self.s0gcm = mesh.project3d(md, vector = self.s0gcm, type = 'node')
            
        return self
    
    def check_consistency(self, md, solution, analyses):
        if 'MasstransportAnalysis' in analyses:
            md = class_utils._check_field(md, fieldname = 'smb.desfac', le = 1, scalar = True)
            md = class_utils._check_field(md, fieldname = 'smb.s0gcm', ge = 0, size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            md = class_utils._check_field(md, fieldname = 'smb.rlaps', ge = 0, scalar = True)
            md = class_utils._check_field(md, fieldname = 'smb.rdl', ge = 0, scalar = True)
            md = class_utils._check_field(md, fieldname = 'smb.dailysnowfall', timeseries = True, allow_nan = False, allow_inf = False, ge = 0)
            md = class_utils._check_field(md, fieldname = 'smb.dailyrainfall', timeseries = True, allow_nan = False, allow_inf = False, ge = 0)
            md = class_utils._check_field(md, fieldname = 'smb.dailydsradiation', timeseries = True, allow_nan = False, allow_inf = False, ge = 0)
            md = class_utils._check_field(md, fieldname = 'smb.dailydlradiation', timeseries = True, allow_nan = False, allow_inf = False, ge = 0)
            md = class_utils._check_field(md, fieldname = 'smb.dailywindspeed', timeseries = True, allow_nan = False, allow_inf = False, ge = 0)
            md = class_utils._check_field(md, fieldname = 'smb.dailypressure', timeseries = True, allow_nan = False, allow_inf = False, ge = 0)
            md = class_utils._check_field(md, fieldname = 'smb.dailyairdensity', timeseries = True, allow_nan = False, allow_inf = False, ge = 0)
            md = class_utils._check_field(md, fieldname='smb.dailyairhumidity', timeseries=True, allow_nan=False, allow_inf=False, ge=0)
            md = class_utils._check_field(md, fieldname='smb.dailytemperature', timeseries=True, allow_nan=False, allow_inf=False, ge=0)
            md = class_utils._check_field(md, fieldname = 'smb.ismethod', scalar = True, values = [0, 1])
            md = class_utils._check_field(md, fieldname = 'smb.isdesertification', scalar = True, values = [0, 1])
            md = class_utils._check_field(md, fieldname = 'smb.isLWDcorrect', scalar = True, values = [0, 1])

            if self.ismethod: # transient mode
                md = class_utils._check_field(md, fieldname = 'smb.desfacElevation', scalar = True, ge = 0)
                md = class_utils._check_field(md, fieldname = 'smb.albedo_scheme', scalar = True, values = [0, 1, 2, 3, 4])
                md = class_utils._check_field(md, fieldname = 'smb.alb_smax', scalar = True, ge = 0)
                md = class_utils._check_field(md, fieldname = 'smb.mask', size = (md.mesh.numberofvertices, ), values = [0, 1, 2])
                md = class_utils._check_field(md, fieldname = 'smb.albedo', allow_nan = False, allow_inf = False, size = (md.mesh.numberofvertices, ))
                md = class_utils._check_field(md, fieldname = 'smb.albedo_snow', allow_nan = False, allow_inf = False, size = (md.mesh.numberofvertices, ))
                md = class_utils._check_field(md, fieldname = 'smb.alb_smax', ge = 0, le = 1, allow_nan = False, allow_inf = False, scalar = True)
                md = class_utils._check_field(md, fieldname = 'smb.alb_smin', ge = 0, le = 1, allow_nan = False, allow_inf = False, scalar = True)
                md = class_utils._check_field(md, fieldname = 'smb.albi', ge = 0, le = 1, allow_nan = False, allow_inf = False, scalar = True)
                md = class_utils._check_field(md, fieldname = 'smb.albl', ge = 0, le = 1, allow_nan = False, allow_inf = False, scalar = True)
                md = class_utils._check_field(md, fieldname = 'smb.hice', allow_nan=False, allow_inf=False, size=(md.mesh.numberofvertices, ))
                md = class_utils._check_field(md, fieldname='smb.hsnow', allow_nan=False, allow_inf=False, size=(md.mesh.numberofvertices, ))
                md = class_utils._check_field(md, fieldname = 'smb.qmr', allow_nan=False, allow_inf=False, size=(md.mesh.numberofvertices, ))
        md = class_utils._check_field(md, fieldname = 'smb.steps_per_step', ge = 1, scalar = True)
        md = class_utils._check_field(md, fieldname = 'smb.averaging', scalar = True, values = [0, 1, 2])
        md = class_utils._check_field(md, fieldname = 'smb.requested_outputs', string_list = True)

        return md
    
    # Initialise empty fields of correct dimensions
    def _initialize(self, md):
        """
        Initialise empty fields in smb.semic.
        """

        if np.isnan(self.s0gcm):
                if hasattr(md.geometry, 'surface') and md.geometry.surface.size == md.mesh.numberofvertices:
                    self.s0gcm = md.geometry.surface
                    print('      no SMBsemic.s0gcm specified: values from md.geometry.surface')
                else:
                    self.s0gcm = np.zeros((md.mesh.numberofvertices,))
                    print('      no SMBsemic.s0gcm specified: values set as zero')

        if np.isnan(self.mask):
            self.mask = 2 * np.ones((md.mesh.numberofvertices, ))
            print('      no SMBsemic.mask specified: values set as 2 for ice')

        self.Tamp = 3 * np.ones((md.mesh.numberofvertices, ))
        self.hice = 10 * np.ones((md.mesh.numberofvertices, ))
        self.hsnow = 5 * np.ones((md.mesh.numberofvertices, ))
        self.qmr = np.zeros((md.mesh.numberofvertices, ))
        self.dailywindspeed = np.zeros((md.mesh.numberofvertices, ))

        return self

    # Process requested outputs, expanding 'default' to appropriate outputs
    def _process_outputs(self,
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
        
        # Marshall method for saving the smb.semic parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [smb.semic] parameters to a binary file.

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
        execute._write_model_field(fid, prefix, name = 'md.smb.model', data = 12, format = 'Integer')

        ## Write DoubleMat fields
        fieldnames = ['dailysnowfall', 'dailyrainfall', 'dailydsradiation', 'dailydlradiation',
                      'dailywindspeed', 'dailypressure', 'dailyairdensity', 'dailyairhumidity', 'dailytemperature']
        for field in fieldnames:
            execute._write_model_field(fid, prefix, obj = self, fieldname = field, format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

        ## Write Integer fields
        fieldnames = ['isdesertification', 'isLWDcorrect', 'steps_per_step', 'averaging']
        for field in fieldnames:
            execute._write_model_field(fid, prefix, obj = self, fieldname = field, format = 'Integer')

        ## Write conditional fields
        if self.ismethod:
            ### Write DoubleMat fields
            fieldnames = ['Tamp', 'mask', 'hice', 'hsnow', 'qmr', 'albedo', 'albedo_snow']
            for field in fieldnames:
                execute._write_model_field(fid, prefix, obj = self, fieldname = field, format = 'DoubleMat', mattype = 1)

            ### Write Double fields
            fieldnames = ['hcrit', 'rcrit', 'alb_smax', 'alb_smin', 'albi', 'albl', 'tmin', 'tmax', 'mcrit', 'wcrit', 'tau_a', 'tau_f', 'tmid', 'afac']
            for field in fieldnames:
                execute._write_model_field(fid, prefix, obj = self, fieldname = field, format = 'Double')
            
            ### Write Integer fields
            execute._write_model_field(fid, prefix, obj = self, fieldname = 'albedo_scheme', format = 'Integer')

        ## Write other fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'ismethod', format = 'Integer')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'desfac', format = 'Double')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'desfacElevation', format = 'Double')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 's0gcm', format = 'DoubleMat', mattype = 1)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'rlaps', format = 'Double')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'rdl', format = 'Double')
        execute._write_model_field(fid, prefix, name = 'md.smb.requested_outputs', data = self._process_outputs(md), format = 'StringArray')