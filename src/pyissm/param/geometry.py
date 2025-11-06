import numpy as np
from . import param_utils
from . import class_registry
from .. import execute

@class_registry.register_class
class geometry(class_registry.manage_state):
    """
    Geometry parameters class for ISSM.

    This class encapsulates geometric parameters that define the ice sheet geometry in the ISSM (Ice Sheet System Model) framework.
    It stores elevation data for ice surface, thickness, base, and bed that are fundamental to ice sheet modeling.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    surface : ndarray, default=nan
        Ice upper surface elevation [m].
    thickness : ndarray, default=nan
        Ice thickness [m].
    base : ndarray, default=nan
        Ice base elevation [m].
    bed : ndarray, default=nan
        Bed elevation [m].
    hydrostatic_ratio : float, default=nan
        Hydrostatic ratio for floating ice.

    Methods
    -------
    __init__(self, other=None)
        Initializes the geometry parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the geometry parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.geometry = pyissm.param.geometry()
    md.geometry.surface = surface_elevation
    md.geometry.thickness = ice_thickness
    md.geometry.bed = bed_elevation
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.surface = np.nan
        self.thickness = np.nan
        self.base = np.nan
        self.bed = np.nan
        self.hydrostatic_ratio = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   geometry parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'surface', 'ice upper surface elevation [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'thickness', 'ice thickness [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'base', 'ice base elevation [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'bed', 'bed elevation [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'hydrostatic_ratio', 'hydrostatic ratio for floating ice'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - geometry Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):  # {{{
        # Early return if LoveSolution
        if solution == 'LoveSolution':
            return md
        else:
            param_utils.check_field(md, fieldname = 'geometry.surface', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            param_utils.check_field(md, fieldname = 'geometry.base', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            param_utils.check_field(md, fieldname = 'geometry.thickness', ge = 0, size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            if any(abs(self.thickness - self.surface + self.base) > 1e-9):
                md.checkmessage('equality thickness = surface-base violated')
            if solution == 'TransientSolution' and md.transient.isgroundingline:
                param_utils.check_field(md, fieldname = 'geometry.bed', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
                if np.any(self.bed - self.base > 1e-12):
                    md.checkmessage('base < bed on one or more vertices')
                pos = np.where(md.mask.ocean_levelset > 0)
                if np.any(np.abs(self.bed[pos] - self.base[pos]) > 1e-9):
                    md.checkmessage('equality base = bed on grounded ice violated')
                param_utils.check_field(md, fieldname = 'geometry.bed', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)

        return md

    # Marshall method for saving the geometry parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [geometry] parameters to a binary file.

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

        ## 1. Handle thickness field
        # Determine the length of the thickness array (could be list or ndarray)
        if isinstance(self.thickness, (list, np.ndarray)):
            length_thickness = len(self.thickness)
        else:
            length_thickness = 1

        # Write thickness data depending on whether it matches number of vertices or elements
        if (length_thickness == md.mesh.numberofvertices) or (length_thickness == md.mesh.numberofvertices + 1):
            execute.WriteData(fid, prefix, obj = self, fieldname = 'thickness', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        elif (length_thickness == md.mesh.numberofelements) or (length_thickness == md.mesh.numberofelements + 1):
            execute.WriteData(fid, prefix, obj = self, fieldname = 'thickness', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofelements + 1, yts = md.constants.yts)
        else:
            # Raise error if thickness does not match expected sizes
            raise RuntimeError('geometry thickness time series should be a vertex or element time series')

        ## 2. Write other geometry fields to file (all fields are of the same type/format)
        fieldnames = list(self.__dict__.keys())
        fieldnames.remove('thickness') # Remove thickness as it is handled separately above
        for fieldname in fieldnames:
                execute.WriteData(fid, prefix, obj = self, fieldname = fieldname, format = 'DoubleMat', mattype = 1)