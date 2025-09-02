import numpy as np
from . import param_utils
from . import class_registry
from .. import execute

@class_registry.register_class
class stochasticforcing(class_registry.manage_state):
    """
    Stochastic forcing parameters class for ISSM.

    This class encapsulates parameters for stochastic forcing in the ISSM (Ice Sheet System Model) framework.
    It allows users to apply random forcing to various physical processes such as surface mass balance,
    basal melting, and calving, enabling uncertainty quantification and probabilistic modeling.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    isstochasticforcing : int, default=0
        Is stochasticity activated?
    fields : str, default='List of fields'
        Fields with stochasticity applied, e.g. ['SMBautoregression'], or ['SMBforcing','DefaultCalving'].
    defaultdimension : int, default=0
        Dimensionality of the noise terms (does not apply to fields with their specific dimension).
    default_id : ndarray, default=nan
        ID of each element for partitioning of the noise terms (does not apply to fields with their specific partition).
    covariance : float, default=nan
        Covariance matrix for within- and between-fields covariance (units must be squared field units), multiple matrices can be concatenated along 3rd dimension to apply different covariances in time.
    timecovariance : float, default=nan
        Starting dates at which covariances apply (only applicable if multiple covariance matrices are prescribed).
    stochastictimestep : float, default=0
        Timestep at which new stochastic noise terms are generated (default: md.timestepping.time_step).
    randomflag : int, default=1
        Whether to apply real randomness (true) or pseudo-randomness with fixed seed (false).

    Methods
    -------
    __init__(self, other=None)
        Initializes the stochasticforcing parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the stochasticforcing parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.stochasticforcing = pyissm.param.stochasticforcing()
    md.stochasticforcing.isstochasticforcing = 1
    md.stochasticforcing.fields = ['SMBforcing']
    md.stochasticforcing.defaultdimension = 1
    md.stochasticforcing.covariance = covariance_matrix
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.isstochasticforcing = 0
        self.fields = 'List of fields'
        self.defaultdimension = 0
        self.default_id = np.nan
        self.covariance = np.nan
        self.timecovariance = np.nan
        self.stochastictimestep = 0
        self.randomflag = 1

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   stochasticforcing parameters:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isstochasticforcing', 'is stochasticity activated?'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'fields', 'fields with stochasticity applied, ex: [\'SMBautoregression\'], or [\'SMBforcing\',\'DefaultCalving\']'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'defaultdimension', 'dimensionality of the noise terms (does not apply to fields with their specific dimension)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'default_id', 'id of each element for partitioning of the noise terms (does not apply to fields with their specific partition)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'covariance', 'covariance matrix for within- and between-fields covariance (units must be squared field units),multiple matrices can be concatenated along 3rd dimension to apply different covariances in time'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'timecovariance', 'starting dates at which covariances apply (only applicabe if multiple covariance matrices are prescribed)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'stochastictimestep', 'timestep at which new stochastic noise terms are generated (default: md.timestepping.time_step)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'randomflag', 'whether to apply real randomness (true) or pseudo-randomness with fixed seed (false)'))
        s += 'Available fields:\n'
        s += '   BasalforcingsDeepwaterMeltingRatearma\n'
        s += '   BasalforcingsSpatialDeepwaterMeltingRate\n'
        s += '   DefaultCalving\n'
        s += '   FloatingMeltRate\n'
        s += '   FrictionWaterPressure\n'
        s += '   FrictionCoulombWaterPressure\n'
        s += '   FrictionSchoofWaterPressure\n'
        s += '   FrontalForcingsRignotarma (thermal forcing)\n'
        s += '   FrontalForcingsSubglacialDischargearma\n'
        s += '   SMBarma\n'
        s += '   SMBforcing\n'
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - stochasticforcing Class'
        return s

    # Marshall method for saving the stochasticforcing parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [stochasticforcing] parameters to a binary file.

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

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'isstochasticforcing', format = 'Boolean')

        ## Conditional writing
        ## NOTE: Taken from $ISSM_DIR/src/classes/stochasticforcing.py
        if self.isstochasticforcing:
            
            ## Initalise
            num_fields = len(self.fields)

            ## Set default timestep
            if self.stochastictimestep == 0:
                self.stochastictimestep = md.timestepping.time_step
            
            if ((type(md.hydrology).__name__ == 'armapw') and md.transient.ishydrology == 1):
                ispwHydroarma = 1
            else:
                ispwHydroarma = 0
            
            ## Get dimensionality of each field
            dimensions = self.defaultdimension * np.ones((num_fields))
            for ind, field in enumerate(self.fields):
                
                ## Check for specific dimensions
                if field == 'SMBarma':
                    dimensions[ind] = md.smb.num_basins
                elif field == 'FrontalForcingsRignotarma':
                    dimensions[ind] = md.frontalforcings.num_basins
                elif field == 'FrontalForcingsSubglacialDischargearma':
                    dimensions[ind] = md.frontalforcings.num_basins
                elif field == 'BasalforcingsDeepwaterMeltingRatearma':
                    dimensions[ind] = md.basalforcings.num_basins
                elif field == 'FrictionWaterPressure' and ispwHydroarma:
                    dimensions[ind] = md.hydrology.num_basins

            if(len(np.shape(self.covariance)) == 3):
                nrow, ncol, numtcovmat = np.shape(self.covariance)
                lsCovmats = []
                for ii in range(numtcovmat):
                    lsCovmats.append(self.covariance[:, :, ii])
                if(md.timestepping.interp_forcing == 1):
                    print('WARNING: md.timestepping.interp_forcing is 1, but be aware that there is no interpolation between covariance matrices')
                    print('         the changes between covariance matrices occur at the time steps specified in md.stochasticforcing.timecovariance')
            elif(len(np.shape(self.covariance)) == 2):
                nrow, ncol = np.shape(self.covariance)
                numtcovmat = 1
                lsCovmats = [self.covariance]
            
            ## Scaling covariance matrix (scale column-by-column and row-by-row)
                ## list of fields that need scaling * 1 / yts
            scaledfields = ['BasalforcingsDeepwaterMeltingRatearma','BasalforcingsSpatialDeepwaterMeltingRate','DefaultCalving', 'FloatingMeltRate', 'SMBarma', 'SMBforcing']
            tempcovariance2d = np.zeros((numtcovmat,nrow*ncol))
            
            ## Loop over covariance matrices
            for kk in range(numtcovmat):
                kkcov = lsCovmats[kk]
                ## Loop over the fields
                for i in range(num_fields):
                    if self.fields[i] in scaledfields:
                        inds = range(int(np.sum(dimensions[0:i])), int(np.sum(dimensions[0:i + 1])))
                        for row in inds:  # scale rows corresponding to scaled field
                            kkcov[row, :] = 1 / md.constants.yts * kkcov[row, :]
                        for col in inds:  # scale columns corresponding to scaled field
                            kkcov[:, col] = 1 / md.constants.yts * kkcov[:, col]
                ## Save scaled covariance
                for rr in range(nrow):
                    ind0 = rr*ncol
                    tempcovariance2d[kk,ind0:ind0+ncol] = np.copy(kkcov[rr,:])
            
            ## Set dummy default_id vector if defaults not used
            if np.any(np.isnan(self.default_id)):
                self.default_id = np.zeros(md.mesh.numberofelements)

            ## Set dummy timecovariance vector if a single covariance matrix is used
            if(numtcovmat==1):
                self.timecovariance = np.array([md.timestepping.start_time])

            ## Reshape dimensions as column array for marshalling
            dimensions = dimensions.reshape(1, len(dimensions))

            ## Write fields
            execute.WriteData(fid, prefix, name = 'md.stochasticforcing.num_fields', data = num_fields, format = 'Integer')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'fields', format = 'StringArray')
            execute.WriteData(fid, prefix, name = 'md.stochasticforcing.dimensions', data = dimensions, format = 'IntMat', mattype = 2)
            execute.WriteData(fid, prefix, name = 'md.stochasticforcing.default_id', data = self.default_id - 1, format = 'IntMat', mattype = 2) # 0-indexed
            execute.WriteData(fid, prefix, obj = self, fieldname = 'defaultdimension', format = 'Integer')
            execute.WriteData(fid, prefix, name = 'md.stochasticforcing.num_timescovariance', data = numtcovmat, format = 'Integer')
            execute.WriteData(fid, prefix, name = 'md.stochasticforcing.covariance', data = tempcovariance2d, format = 'DoubleMat')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'timecovariance', format = 'DoubleMat', scale = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'stochastictimestep', format = 'Double', scale = md.constants.yts)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'randomflag', format = 'Boolean')