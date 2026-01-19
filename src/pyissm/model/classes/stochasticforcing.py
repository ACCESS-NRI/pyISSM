import numpy as np
import warnings

from pyissm.model.classes import friction
from pyissm.model.classes import hydrology
from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm.model import execute, mesh

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
    md.stochasticforcing = pyissm.model.classes.stochasticforcing()
    md.stochasticforcing.isstochasticforcing = 1
    md.stochasticforcing.fields = ['SMBforcing']
    md.stochasticforcing.defaultdimension = 1
    md.stochasticforcing.covariance = covariance_matrix
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.isstochasticforcing = 0
        self.fields = []
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
        s += '{}\n'.format(class_utils.fielddisplay(self, 'isstochasticforcing', 'is stochasticity activated?'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'fields', 'fields with stochasticity applied, ex: [\'smb.default\',\'calving.default\']'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'defaultdimension', 'dimensionality of the noise terms (does not apply to fields with their specific dimension)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'default_id', 'id of each element for partitioning of the noise terms (does not apply to fields with their specific partition)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'covariance', 'covariance matrix for within- and between-fields covariance (units must be squared field units),multiple matrices can be concatenated along 3rd dimension to apply different covariances in time'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'timecovariance', 'starting dates at which covariances apply (only applicabe if multiple covariance matrices are prescribed)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'stochastictimestep', 'timestep at which new stochastic noise terms are generated (default: md.timestepping.time_step)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'randomflag', 'whether to apply real randomness (true) or pseudo-randomness with fixed seed (false)'))
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
    
    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude stochasticforcing fields to 3D
        """
        self.default_id = mesh.project_3d(md, vector = self.default_id, type = 'element')

        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        # Early return if stochasticforcing not enabled
        if not self.isstochasticforcing:
            return md

        num_fields = len(self.fields)

        # If not specified, set stochastictimestep to md.timestepping.time_step
        if self.stochastictimestep == 0:
            md.stochasticforcing.stochastictimestep = md.timestepping.time_step
            warnings.warn('pyissm.model.classes.stochasticforcing.check_consistency: stochasticforcing.stochastictimestep not specified -- set to md.timestepping.time_step')

        # Covariance matrix checks
        if(len(np.shape(self.covariance))==3):
           # Set number of covariance matrices in time
           numtcovmat = np.shape(self.covariance)[2]
           lsCovmats = []
           for ii in range(numtcovmat):
               lsCovmats.append(self.covariance[:,:,ii])
               try:
                   np.linalg.cholesky(self.covariance[:,:,ii])
               except:
                   raise TypeError('pyissm.model.classes.stochasticforcing.check_consistency: an entry in md.stochasticforcing.covariance is not positive definite')
        elif(len(np.shape(self.covariance))==2):
            numtcovmat = 1
            lsCovmats = [self.covariance]
            # Check that covariance matrix is positive definite (this is done internally by linalg)
            try:
                np.linalg.cholesky(self.covariance)
            except:
                raise TypeError('pyissm.model.classes.stochasticforcing.check_consistency:md.stochasticforcing.covariance is not positive definite')

        # Check that all fields agree with the corresponding md class and if any field needs the default params
        checkdefaults = False  # Need to check defaults only if one of the fields does not have its own dimensionality
        structstoch = class_utils.supported_stochastic_forcings(return_dict = True)
        
        # Check if hydrologyarmapw is used
        if isinstance(md.hydrology, hydrology.armapw) and md.transient.ishydrology==1:
            ispwHydroarma = 1
        else:
            ispwHydroarma = 0

        component_map = {
            'SMB': md.smb,
            'FrontalForcings': md.frontalforcings,
            'Calving': md.calving,
            'Basalforcings': md.basalforcings,
            'Friction': md.friction
        }

        for field in self.fields:
            expected_class = structstoch.get(field)
            if expected_class is None:
                raise ValueError(f'Field {field} in stochasticforcing is not recognized.')
            
            # Identify which md submodel this field belongs to
            component = next((obj for key, obj in component_map.items() if key in field), None)
            if component is None:
                raise ValueError(f"Cannot identify model component for stochasticforcing field '{field}'")
            
            # Handle friction (WaterPressure)
            if "WaterPressure" in field:
                if not isinstance(md.friction, (friction.default, friction.coulomb, friction.schoof)):
                    raise TypeError(f"md.friction does not agree with stochasticforcing field '{field}'")

                if md.friction.coupling not in [0, 1, 2]:
                    raise TypeError(
                        f"stochasticforcing field '{field}' only supported for md.friction.coupling in [0, 1, 2]"
                    )

                if isinstance(md.friction, friction) and np.any(md.friction.q == 0):
                    raise TypeError(f"stochasticforcing field '{field}' requires non-zero q exponent")

                continue

            # Handle other components
            if not isinstance(component, expected_class):
                raise TypeError(
                    f"{component.__class__.__name__} does not agree with stochasticforcing field '{field}' "
                    f"(expected instance of {expected_class.__name__})"
                    )

            # Check dimensions
            if (field not in ["SMBarma", "FrontalForcingsRignotarma", "BasalforcingsDeepwaterMeltingRatearma"] and not (field == "FrictionWaterPressure" and ispwHydroarma)):
                checkdefaults = True 


        # Retrieve sum of all the field dimensionalities
        dimensions = self.defaultdimension * np.ones((num_fields))

        # Mapping of supported ARMA fields to their md components
        # NOTE: Maintain legacy naming for compatibility with MATLAB version
        arma_fields = {
            "SMBarma": ("smb", "arma_timestep", "num_basins", "indSMBarma"),
            "FrontalForcingsRignotarma": ("frontalforcings", "arma_timestep", "num_basins", "indTFarma"),
            "FrontalForcingsSubglacialDischargearma": ("frontalforcings", "sd_arma_timestep", "num_basins", "indSdarma"),
            "BasalforcingsDeepwaterMeltingRatearma": ("basalforcings", "arma_timestep", "num_basins", "indBDWarma"),
            "hydrologyarmapw": ("hydrology", "arma_timestep", "num_basins", "indPwarma"),
        }

        # Initialise index variables for compatibility
        indSMBarma = indTFarma = indSdarma = indBDWarma = indPwarma = -1
        indices = {}
        timesteps = {}

        # Assign dimensions and check timesteps
        for field, (attr, ts_attr, dim_attr) in arma_fields.items():
            if field in self.fields:
                idx = self.fields.index(field)
                obj = getattr(md, attr)
                dimensions[idx] = getattr(obj, dim_attr)
                indices[field] = idx
                timesteps[field] = getattr(obj, ts_attr)
                if timesteps[field] < self.stochastictimestep:
                    raise TypeError(f"{field} cannot have a timestep shorter than stochastictimestep")

        size_tot = np.sum(dimensions)

        # Check pairwise covariance consistency
        present_fields = list(indices.keys())
        for i, f1 in enumerate(present_fields):
            for f2 in present_fields[i + 1 :]:
                ts1, ts2 = timesteps[f1], timesteps[f2]
                if ts1 == ts2:
                    continue  # same timestep, OK

                i1, i2 = indices[f1], indices[f2]
                for covm in lsCovmats:
                    r1, r2 = int(np.sum(dimensions[:i1])), int(np.sum(dimensions[:i1 + 1]))
                    c1, c2 = int(np.sum(dimensions[:i2])), int(np.sum(dimensions[:i2 + 1]))
                    covsum = covm[r1:r2, c1:c2]
                    if np.any(covsum != 0):
                        raise IOError(
                            f"{f1} and {f2} have different arma_timestep and non-zero covariance"
                        )

        class_utils.check_field(md, fieldname = 'stochasticforcing.isstochasticforcing', values = [0, 1])
        class_utils.check_field(md, fieldname = 'stochasticforcing.fields', numel = num_fields, cell = True, values = class_utils.supported_stochastic_forcings())
        class_utils.check_field(md, fieldname = 'stochasticforcing.covariance', size = (size_tot, size_tot, numtcovmat), allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'stochasticforcing.stochastictimestep', ge = md.timestepping.time_step, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'stochasticforcing.randomflag', scalar = True, values = [0, 1])
        if(numtcovmat > 1):
            class_utils.check_field(md, fieldname = 'stochasticforcing.timecovariance', ge = md.timestepping.start_time, le = md.timestepping.final_time, size = (1, numtcovmat), allow_nan = False, allow_inf = False)
        if checkdefaults:
            class_utils.check_field(md, fieldname = 'stochasticforcing.defaultdimension', scalar = True, gt = 0, allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = 'stochasticforcing.default_id', ge = 0, le = self.defaultdimension, size = (md.mesh.numberofelements, ), allow_nan = False, allow_inf = False)

        return md

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