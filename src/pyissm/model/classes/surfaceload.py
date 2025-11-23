import numpy as np
from pyissm.param import param_utils
from pyissm.param import class_registry
from pyissm import execute

@class_registry.register_class
class surfaceload(class_registry.manage_state):
    """
    Surface load parameters class for ISSM.

    This class encapsulates parameters for surface loading in the ISSM (Ice Sheet System Model) framework.
    Surface loads include ice thickness changes, water height changes, and other surface loads (e.g., sediments)
    that affect solid Earth deformation and sea level calculations.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    icethicknesschange : ndarray, default=nan
        Thickness change: ice height equivalent [mIce/yr].
    waterheightchange : ndarray, default=nan
        Water height change: water height equivalent [mWater/yr].
    other : ndarray, default=nan
        Other loads (sediments) [kg/m^2/yr].

    Methods
    -------
    __init__(self, other=None)
        Initializes the surfaceload parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the surfaceload parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.surfaceload = pyissm.param.surfaceload()
    md.surfaceload.icethicknesschange = ice_thickness_change
    md.surfaceload.waterheightchange = water_height_change
    md.surfaceload.other = sediment_load
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.icethicknesschange = np.nan
        self.waterheightchange = np.nan
        self.otherchange = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   surfaceload:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'icethicknesschange', 'thickness change: ice height equivalent [mIce/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'waterheightchange', 'water height change: water height equivalent [mWater/yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'otherchange', 'other loads (sediments) [kg/m^2/yr]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - surfaceload Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):  # {{{
        # Early return if required analysis/solutions are not present
        if ('SealevelchangeAnalysis' not in analyses) or (solution == 'TransientSolution' and not md.transient.isslc):
            return md
        
        if type(self.icethicknesschange) == np.ndarray:
            param_utils.check_field(md, fieldname = 'solidearth.surfaceload.icethicknesschange', timeseries = True, allow_nan = False, allow_inf = False)
        if type(self.waterheightchange) == np.ndarray:
            param_utils.check_field(md, fieldname = 'solidearth.surfaceload.waterheightchange', timeseries = True, allow_nan = False, allow_inf = False)
        if type(self.otherchange) == np.ndarray:
            param_utils.check_field(md, fieldname = 'solidearth.surfaceload.otherchange', timeseries = True, allow_nan = False, allow_inf = False)

        return md
    
    # Marshall method for saving the surfaceload parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [surfaceload] parameters to a binary file.

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

        ## Check fields
        if np.isnan(self.icethicknesschange):
            self.icethicknesschange = np.zeros((md.mesh.numberofelements + 1, ))

        if np.isnan(self.waterheightchange):
            self.waterheightchange = np.zeros((md.mesh.numberofelements + 1, ))

        if np.isnan(self.otherchange):
            self.otherchange = np.zeros((md.mesh.numberofelements + 1, ))

        ## Write fields
        execute.WriteData(fid, prefix, name = 'md.solidearth.surfaceload.icethicknesschange', data = self.icethicknesschange, format = 'MatArray', timeserieslength = md.mesh.numberofelements + 1, yts = md.constants.yts, scale = 1 / md.constants.yts)
        execute.WriteData(fid, prefix, name = 'md.solidearth.surfaceload.waterheightchange', data = self.waterheightchange, format = 'MatArray', timeserieslength = md.mesh.numberofelements + 1, yts = md.constants.yts, scale = 1 / md.constants.yts)
        execute.WriteData(fid, prefix, name = 'md.solidearth.surfaceload.otherchange', data = self.otherchange, format = 'MatArray', timeserieslength = md.mesh.numberofelements + 1, yts = md.constants.yts, scale = 1 / md.constants.yts)

