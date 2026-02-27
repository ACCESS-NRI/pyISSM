"""
Dynamic Sea Level (DSL) classes for ISSM.
"""

import numpy as np
import warnings
from pyissm.model.classes import class_utils, class_registry
from pyissm.model import execute, mesh

## ------------------------------------------------------
## dsl.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    """
    Default DSL class for ISSM.

    This class contains the default parameters for dynamic sea level (DSL) in the ISSM framework.
    It defines the main DSL-related parameters.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    global_average_thermosteric_sea_level : :class:`float` or :class:`numpy.ndarray`, default=np.nan
        Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m).
    sea_surface_height_above_geoid : :class:`float` or :class:`numpy.ndarray`, default=np.nan
        Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m).
    sea_water_pressure_at_sea_floor : :class:`float` or :class:`numpy.ndarray`, default=np.nan
        Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!).

    Examples
    --------
    .. code-block:: python
    
        >>> md.dsl = pyissm.model.classes.dsl.default()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.global_average_thermosteric_sea_level = np.nan
        self.sea_surface_height_above_geoid = np.nan
        self.sea_water_pressure_at_sea_floor = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   dsl parameters:\n'
        s += '{}\n'.format(class_utils.fielddisplay(self, 'global_average_thermosteric_sea_level', 'Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m).'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'sea_surface_height_above_geoid', 'Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m).'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'sea_water_pressure_at_sea_floor', 'Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!).'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - dsl Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude dsl.default fields to 3D
        """
        self.sea_surface_height_above_geoid = mesh.project_3d(md, vector = self.sea_surface_height_above_geoid, type = 'node', layer = 1)
        self.sea_water_pressure_at_sea_floor = mesh.project_3d(md, vector = self.sea_water_pressure_at_sea_floor, type = 'node', layer = 1)

        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [dsl.default] parameters.

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

        # Early return if no sealevelchange analysis or if transient solution without isslc or oceantransport
        if ('SealevelchangeAnalysis' not in analyses) or (solution == 'TransientSolution' and not md.transient.isslc) or (not md.transient.isoceantransport):
            return md
        class_utils.check_field(md, fieldname = 'dsl.global_average_thermosteric_sea_level', allow_nan = True, allow_inf = True)
        class_utils.check_field(md, fieldname = 'dsl.sea_surface_height_above_geoid', allow_nan = True, allow_inf = True, timeseries = True)
        class_utils.check_field(md, fieldname = 'dsl.sea_water_pressure_at_sea_floor', allow_nan = True, allow_inf = True, timeseries = True)

        if md.solidearth.settings.compute_bp_grd:
            class_utils.check_field(md, fieldname = 'dsl.sea_water_pressure_at_sea_floor', allow_empty = True)

        return md
    
    # Initialise empty fields of correct dimensions
    def initialise(self, md):
        """
        Initialise [dsl.default] empty fields.

        If current values of required fields are np.nan, they will be set to default required shapes/values and warnings will be issued.

        Examples
        --------
        .. code-block:: python

            >>> md.dsl = pyissm.model.classes.dsl.default()
            # At this point, initial fields are np.nan
            # After calling initialise, they will be set to default shapes/values with warnings issued.
            >>> md.dsl.initialise(md)
        """

        if np.all(np.isnan(self.global_average_thermosteric_sea_level)):
            self.global_average_thermosteric_sea_level = np.array([0, 0]).reshape(-1, 1)
            warnings.warn('pyissm.model.classes.dsl.default: no dsl.global_average_thermosteric_sea_level specified: transient values set to zero')

        if np.all(np.isnan(self.sea_surface_height_above_geoid)):
            self.sea_surface_height_above_geoid = np.append(np.zeros(md.mesh.numberofvertices, ), 0).reshape(-1, 1)
            warnings.warn('pyissm.model.classes.dsl.default: no dsl.sea_surface_height_above_geoid specified: transient values set to zero')

        if np.all(np.isnan(self.sea_water_pressure_at_sea_floor)):
            self.sea_water_pressure_at_sea_floor = np.append(np.zeros(md.mesh.numberofvertices, ), 0).reshape(-1, 1)
            warnings.warn('pyissm.model.classes.dsl.default: no dsl.sea_water_pressure_at_sea_floor specified: transient values set to zero')

        return self

    # Marshall method for saving the dsl.default parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [dsl.default] parameters to a binary file.

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

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.dsl.model', data = 1, format = 'Integer')

        ## Write DoubleMat fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'global_average_thermosteric_sea_level', format = 'DoubleMat', mattype = 2, timeserieslength = 2, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'sea_surface_height_above_geoid', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'sea_water_pressure_at_sea_floor', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)


## ------------------------------------------------------
## dsl.mme
## ------------------------------------------------------
@class_registry.register_class
class mme(class_registry.manage_state):
    """
    Multi-Model Ensemble (MME) DSL class for ISSM.

    This class conatins the parameters for dynamic sea level (DSL) for a multi-model ensemble (MME) of CMIP5 outputs in the ISSM framework.
    It defines the main DSL-related parameters for each ensemble member.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    modelid : :class:`int`, default=0
        Index into the multi-model ensemble, determines which field will be used.
    global_average_thermosteric_sea_level : :class:`float` or :class:`numpy.ndarray`, default=np.nan
        Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m) for each ensemble member.
    sea_surface_height_above_geoid : :class:`float` or :class:`numpy.ndarray`, default=np.nan
        Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m) for each ensemble member.
    sea_water_pressure_at_sea_floor : :class:`float` or :class:`numpy.ndarray`, default=np.nan
        Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!) for each ensemble member.

    Examples
    --------
    .. code-block:: python
        
        >>> md.dsl = pyissm.model.classes.dsl.mme()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.modelid = 0
        self.global_average_thermosteric_sea_level = []
        self.sea_surface_height_above_geoid = []
        self.sea_water_pressure_at_sea_floor = []

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   dsl mme parameters:\n'
        s += '{}\n'.format(class_utils.fielddisplay(self, 'modelid', 'index into the multi-model ensemble, determines which field will be used.'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'global_average_thermosteric_sea_level', 'Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m) for each ensemble.'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'sea_surface_height_above_geoid', 'Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m) for each ensemble.'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'sea_water_pressure_at_sea_floor', 'Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!) for each ensemble.'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - dsl mme Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude dsl.mme fields to 3D
        """
        for i in range(len(self.global_average_thermosteric_sea_level)):
            self.sea_surface_height_above_geoid[i] = mesh.project_3d(md, vector = self.self.sea_surface_height_above_geoid[i], type = 'node', layer = 1)
            self.sea_water_pressure_at_sea_floor[i] = mesh.project_3d(md, vector = self.sea_water_pressure_at_sea_floor[i], type = 'node', layer = 1)

        return self    

    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [dsl.mme] parameters.

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

        # Early return if no sealevelchange analysis or if transient solution without isslc or oceantransport
        if ('SealevelchangeAnalysis' not in analyses) or (solution == 'TransientSolution' and not md.transient.isslc) or (not md.transient.isoceantransport):
            return md

        for i in range(len(self.global_average_thermosteric_sea_level)):
            class_utils.check_field(md, field = self.global_average_thermosteric_sea_level[i], allow_nan = True, allow_inf = True)
            class_utils.check_field(md, field = self.sea_surface_height_above_geoid[i], allow_nan = True, allow_inf = True, timeseries = True)
            class_utils.check_field(md, field = self.sea_water_pressure_at_sea_floor[i], allow_nan = True, allow_inf = True, timeseries = True)
        
        class_utils.check_field(md, field = self.modelid, allow_nan = True, allow_inf = True, ge = 1, le = len(self.global_average_thermosteric_sea_level))

        if self.solidearth.settings.compute_bp_grd:
            class_utils.check_field(md, field = self.sea_water_pressure_at_sea_floor, allow_empty = True)

        return md
    
    # Marshall method for saving the dsl.default parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [dsl.mme] parameters to a binary file.

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

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.dsl.model', data = 2, format = 'Integer')

        ## Write DoubleMat fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'global_average_thermosteric_sea_level', format = 'MatArray', timeserieslength = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'sea_surface_height_above_geoid', format = 'MatArray', timeserieslength = md.mesh.numberofvertices + 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'sea_water_pressure_at_sea_floor', format = 'MatArray', timeserieslength = md.mesh.numberofvertices + 1)

        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'modelid', format = 'Double')
        execute.WriteData(fid, prefix, name = 'md.dsl.nummodels', data = len(self.global_average_thermosteric_sea_level), format = 'Integer')