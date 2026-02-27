import numpy as np
from pyissm.model.classes import class_utils, class_registry
from pyissm.model import execute, mesh

@class_registry.register_class
class geometry(class_registry.manage_state):
    """
    Geometry class for ISSM.

    This class contains geometric parameters that define the ice sheet geometry in the ISSM framework.
    It stores elevation data for ice surface, thickness, base, and bed.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    surface : :class:`ndarray`, default=np.nan
        Ice upper surface elevation [m].
    thickness : :class:`ndarray`, default=np.nan
        Ice thickness [m].
    base : :class:`ndarray`, default=np.nan
        Ice base elevation [m].
    bed : :class:`numpy.ndarray`, default=np.nan
        Bed elevation [m].
    hydrostatic_ratio : :class:`float`, default=nan
        Hydrostatic ratio for floating ice.

    Examples
    --------
    .. code-block:: python

        >>> md.geometry = pyissm.model.classes.geometry()
        >>> md.geometry.surface = surface_elevation
        >>> md.geometry.thickness = ice_thickness
        >>> md.geometry.bed = bed_elevation
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

        s += '{}\n'.format(class_utils._field_display(self, 'surface', 'ice upper surface elevation [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'thickness', 'ice thickness [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'base', 'ice base elevation [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'bed', 'bed elevation [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'hydrostatic_ratio', 'hydrostatic ratio for floating ice'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - geometry Class'
        return s

    def _extrude(self, md):
        """
        Extrude geometry fields to 3D mesh
        """
        self.surface = mesh._project_3d(md, vector = self.surface, type = 'node')
        self.thickness = mesh._project_3d(md, vector = self.thickness, type = 'node')
        self.hydrostatic_ratio = mesh._project_3d(md, vector = self.hydrostatic_ratio, type = 'node')
        self.base = mesh._project_3d(md, vector = self.base, type = 'node')
        self.bed = mesh._project_3d(md, vector = self.bed, type = 'node')
        
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [geometry] parameters.

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

        # Early return if LoveSolution
        if solution == 'LoveSolution':
            return md
        else:
            class_utils.check_field(md, fieldname = 'geometry.surface', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = 'geometry.base', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = 'geometry.thickness', ge = 0, size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            if any(abs(self.thickness - self.surface + self.base) > 1e-9):
                md.check_message('equality thickness = surface-base violated')
            if solution == 'TransientSolution' and md.transient.isgroundingline:
                class_utils.check_field(md, fieldname = 'geometry.bed', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
                if np.any(self.bed - self.base > 1e-12):
                    md.check_message('base < bed on one or more vertices')
                pos = np.where(md.mask.ocean_levelset > 0)
                if np.any(np.abs(self.bed[pos] - self.base[pos]) > 1e-9):
                    md.check_message('equality base = bed on grounded ice violated')
                class_utils.check_field(md, fieldname = 'geometry.bed', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)

        return md

    # Marshall method for saving the geometry parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [geometry] parameters to a binary file.

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
        fieldnames = ['surface', 'base', 'bed', 'hydrostatic_ratio']
        for fieldname in fieldnames:
                execute.WriteData(fid, prefix, obj = self, fieldname = fieldname, format = 'DoubleMat', mattype = 1)