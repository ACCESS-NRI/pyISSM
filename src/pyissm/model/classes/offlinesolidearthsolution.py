import numpy as np
from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm.model import execute

@class_registry.register_class
class offlinesolidearthsolution(class_registry.manage_state):
    """
    OfflineSolidEarthSolution class for ISSM.

    This class defines the default parameters for the offline solid-Earth solution used in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    displacementeast : float or ndarray, default=np.nan
        Solid-Earth Eastwards bedrock displacement time series (m)
    displacementnorth : float or ndarray, default=np.nan
        Solid-Earth Northwards bedrock displacement time series (m)
    displacementup : float or ndarray, default=np.nan
        Solid-Earth bedrock uplift time series (m)
    geoid : float or ndarray, default=np.nan
        Solid-Earth geoid time series (m)

    Methods
    -------
    __init__(self, other=None)
        Initializes the default parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the object.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.solidearth = pyissm.model.classes.offlinesolidearthsolution()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.displacementeast = np.nan
        self.displacementnorth = np.nan
        self.displacementup = np.nan
        self.geoid = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    def __repr__(self):
        s = '         units for time series is (yr)\n       external: offlinesolidearth solution\n'
        s += '{}\n'.format(class_utils.fielddisplay(self, 'displacementeast', 'solid-Earth Eastwards bedrock displacement time series (m)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'displacementnorth', 'solid-Earth Northwards bedrock displacement time series (m)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'displacementup', 'solid-Earth bedrock uplift time series (m)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'geoid', 'solid-Earth geoid time series (m)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - offlinesolidearthsolution Class'
        return s

    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        # Early return if required analyses and solutions not present
        if ('SealevelchangeAnalysis' not in analyses) or ((solution=='TransientSolution') and (md.solidearth.settings.isgrd==1)): 
            print('pyissm.model.classes.offlinesolidearthsolution.check_consistency: trying to run GRD patterns while supplying an offline solution for those patterns!')
            return md 
        class_utils.check_field(md, fieldname = 'solidearth.external.displacementeast', timeseries = True, allow_inf = False)
        class_utils.check_field(md, fieldname = 'solidearth.external.displacementnorth', timeseries = True, allow_inf = False)
        class_utils.check_field(md, fieldname = 'solidearth.external.displacementup',  timeseries = True, allow_inf = False)
        class_utils.check_field(md, fieldname = 'solidearth.external.geoid', timeseries = True, allow_inf = False)

        return md

    # Marshall method for saving the offlinesolidearthsolution parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [offlinesolidearthsolution] parameters to a binary file.

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

        # Transform time series into time series rates
        # NOTE: Taken from $ISSM_DIR/src/m/classes/offlinesolidearthsolution.py
        if len(np.shape(self.displacementeast)) == 1:
            print('Warning: offlinesolidearthsolution::marshall_class: only one time step provided, assuming the values are rates per year')
            displacementeast_rate = np.append(np.array(self.displacementeast).reshape(-1, 1), 0)
            displacementnorth_rate = np.append(np.array(self.displacementnorth).reshape(-1, 1), 0)
            displacementup_rate = np.append(np.array(self.displacementup).reshape(-1, 1), 0)
            geoid_rate = np.append(np.array(self.geoid).reshape(-1, 1), 0)
        else:
            time = self.displacementeast[-1, :]
            dt = np.diff(time, axis=0)
            displacementeast_rate = np.diff(self.displacementeast[0:-2, :], 1, 1) / dt
            displacementeast_rate = np.append(displacementeast_rate,time[:-1].reshape(1,-1),axis=0)
            displacementnorth_rate = np.diff(self.displacementnorth[0:-2, :], 1, 1) / dt
            displacementnorth_rate = np.append(displacementnorth_rate,time[:-1].reshape(1,-1),axis=0)
            displacementup_rate = np.diff(self.displacementup[0:-2, :], 1, 1) / dt
            displacementup_rate = np.append(displacementup_rate,time[:-1].reshape(1,-1),axis=0)
            geoid_rate = np.diff(self.geoid[0:-2, :], 1, 1) / dt
            geoid_rate = np.append(geoid_rate,time[:-1].reshape(1,-1),axis=0)

        ## Write headers to file
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.solidearth.external.nature', data = 2, format = 'Integer')

        ## Write fields
        execute.WriteData(fid, prefix, name =  'md.solidearth.external.displacementeast', data = displacementeast_rate, format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, name = 'md.solidearth.external.displacementup', data = displacementup_rate, format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, name = 'md.solidearth.external.displacementnorth', data = displacementnorth_rate, format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, name = 'md.solidearth.external.geoid', data = geoid_rate, format = 'DoubleMat',  mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

