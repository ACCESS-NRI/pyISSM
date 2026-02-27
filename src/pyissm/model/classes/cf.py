"""
cf classes for ISSM.
"""

import numpy as np
from pyissm.model.classes import class_utils, class_registry
from pyissm.model import execute, mesh

## ------------------------------------------------------
## cf.levelsetmisfit
## ------------------------------------------------------
@class_registry.register_class
class levelsetmisfit(class_registry.manage_state):
    """
    Level-set misfit response definition class for ISSM.

    This response is commonly used to measure misfit between a modeled level-set
    field (e.g. calving-front / ice-mask level set) and an observed level-set
    field. It supports optional per-vertex weights and an associated data time
    (years from start) indicating when the observation applies.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    name : :class:`str`
        Identifier for this response definition.
    definitionstring : :class:`str`
        Unique output definition identifier, typically ``"OutputdefinitionN"``.
        MATLAB allows ``Outputdefinition1``..``Outputdefinition2000``.
    model_string : :class:`str`
        String key for the modeled field to be compared.
        Example: ``"MaskIceLevelset"``.
    observation : :class:`numpy.ndarray` or :class:`float`
        Observed field to compare against (often in ISSM time-series format).
        Written as a DoubleMat time series with ``timeserieslength = nv + 1``.
        Default is NaN (unset).
    observation_string : :class:`str`
        Identifier/name for the observation dataset.
    weights : :class:`numpy.ndarray` or :class:`float`
        Weight coefficients (per-vertex, often time-series formatted).
        Written as a DoubleMat time series with ``timeserieslength = nv + 1``.
        Default is NaN (unset).
    weights_string : :class:`str`
        Identifier/name for the weights dataset.
    datatime : :class:`float`
        Time in years from start associated with the data. MATLAB writes this as
        ``round(datatime * yts)`` to the binary file.

    Examples
    --------
    .. code-block:: python

        >>> cf = cflevelsetmisfit()
        >>> cf.name = "CalvingFrontPosition"
        >>> cf.definitionstring = "Outputdefinition1"
        >>> cf.model_string = "MaskIceLevelset"
        >>> cf.observation_string = "LevelsetObservations"
        >>> cf.observation = md.mask.ice_levelset
        >>> cf.weights = np.ones(md.mesh.numberofvertices, )
        >>> cf.weights_string = "WeightsLevelsetObservations"
        >>> cf.datatime = time
        >>> md.outputdefinition.definitions.append(cf)
    """
    
    # Define class string for registry and marshalling
    @classmethod
    def issm_enum_string(cls) -> str:
        return "Cflevelsetmisfit"

    # Initialise with default parameters
    def __init__(self, other=None):
        self.name = ""
        self.definitionstring = ""
        self.model_string = ""
        self.observation = np.nan
        self.observation_string = ""
        self.weights = np.nan
        self.weights_string = ""
        self.datatime = 0.0

        # Inherit matching fields from provided class
        super().__init__(other)
    # Define repr
    def __repr__(self):
        s = "   levelsetmisfit:\n"
        s += '{}\n'.format(class_utils._field_display(self, 'name', 'identifier for this levelsetmisfit response'))
        s += '{}\n'.format(class_utils._field_display(self, 'definitionstring', 'unique output definition string (e.g. Outputdefinition1)'))
        s += '{}\n'.format(class_utils._field_display(self, 'model_string', 'string for field that is modeled'))
        s += '{}\n'.format(class_utils._field_display(self, 'observation', 'observed field compared against the model'))
        s += '{}\n'.format(class_utils._field_display(self, 'observation_string', 'identifier for observed field'))
        s += '{}\n'.format(class_utils._field_display(self, 'weights', 'weights applied to the misfit'))
        s += '{}\n'.format(class_utils._field_display(self, 'weights_string', 'identifier for weights'))
        s += '{}\n'.format(class_utils._field_display(self, 'datatime', 'time (years from start) for data-model comparison'))
        return s

    # Define class string
    def __str__(self):
        return "ISSM - cf.levelsetmisfit Class"

    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [cf.levelsetmisfit] fields to 3D
        """

        # weights
        if np.size(self.weights) > 1:
            if not np.all(np.isnan(self.weights)):
                self.weights = mesh._project_3d(md, vector=self.weights, type="node")
        else:
            if not np.isnan(self.weights):
                self.weights = mesh._project_3d(md, vector=self.weights, type="node")

        # observation
        if np.size(self.observation) > 1:
            if not np.all(np.isnan(self.observation)):
                self.observation = mesh._project_3d(md, vector=self.observation, type="node")
        else:
            if not np.isnan(self.observation):
                self.observation = mesh._project_3d(md, vector=self.observation, type="node")

        return self

    # Check model consistency
    def check_consistency(self, md, solution=None, analyses=None):
        """
        Check consistency of the [cf.levelsetmisfit] parameters.

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

        if not isinstance(self.name, str):
            raise TypeError("levelsetmisfit: 'name' must be a string")

        # MATLAB allowed Outputdefinition1..2000
        outputdef_allowed = [f"Outputdefinition{i}" for i in range(1, 2001)]

        class_utils.check_field(
            md,
            fieldname="levelsetmisfit.definitionstring",
            field=self.definitionstring,
            values=outputdef_allowed,
        )

        class_utils.check_field(
            md,
            fieldname="levelsetmisfit.observation",
            field=self.observation,
            timeseries=True,
            allow_nan=True,
            allow_inf=True,
        )
        class_utils.check_field(
            md,
            fieldname="levelsetmisfit.weights",
            field=self.weights,
            timeseries=True,
            allow_nan=True,
            allow_inf=True,
        )
        class_utils.check_field(
            md,
            fieldname="levelsetmisfit.datatime",
            field=self.datatime,
            le=md.timestepping.final_time,
        )

        return md
    
    # Marshall method for saving the cf.levelsetmisfit parameters
    def marshall_class(self, fid, prefix, md=None):
        """
        Marshall [cf.levelsetmisfit] parameters to a binary file.

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

        nv = md.mesh.numberofvertices
        yts = md.constants.yts

        execute._write_model_field(fid, prefix, data=self.name, name="md.cflevelsetmisfit.name", format="String")
        execute._write_model_field(
            fid,
            prefix,
            data=self.definitionstring,
            name="md.cflevelsetmisfit.definitionstring",
            format="String",
        )
        execute._write_model_field(
            fid,
            prefix,
            data=self.model_string,
            name="md.cflevelsetmisfit.model_string",
            format="String",
        )

        execute._write_model_field(
            fid,
            prefix,
            data=self.observation,
            name="md.cflevelsetmisfit.observation",
            format="DoubleMat",
            mattype=1,
            timeserieslength=nv + 1,
            yts=yts,
        )
        execute._write_model_field(
            fid,
            prefix,
            data=self.observation_string,
            name="md.cflevelsetmisfit.observation_string",
            format="String",
        )

        execute._write_model_field(
            fid,
            prefix,
            data=self.weights,
            name="md.cflevelsetmisfit.weights",
            format="DoubleMat",
            mattype=1,
            timeserieslength=nv + 1,
            yts=yts,
        )
        execute._write_model_field(
            fid,
            prefix,
            data=self.weights_string,
            name="md.cflevelsetmisfit.weights_string",
            format="String",
        )

        execute._write_model_field(
            fid,
            prefix,
            data=float(np.round(self.datatime * yts)),
            name="md.cflevelsetmisfit.datatime",
            format="Double",
        )

## ------------------------------------------------------
## cf.surfacesquare
## ------------------------------------------------------
@class_registry.register_class
class surfacesquare(class_registry.manage_state):
    """
    Surface-square cost-function (response) definition class for ISSM.

    Python/pyISSM equivalent of ISSM MATLAB's ``cfsurfacesquare`` class. This
    response is typically used to compare a modeled surface (or base) field
    against observed data with per-vertex weights, optionally as a time series.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    name : :class:`str`
        Identifier for this response definition.
    definitionstring : :class:`str`
        Unique output definition identifier, typically ``"OutputdefinitionN"``.
    surfaceid : :class:`int`
        Which surface to evaluate:
        - 1: surface
        - 2: base
    model_string : :class:`str`
        Name of the modeled field to compare (string key used by ISSM).
    observation : :class:`numpy.ndarray` or :class:`float`
        Observed field to compare against. May be provided as a time series in
        ISSM “timeseries” format (length = numberofvertices + 1).
        Default is NaN (unset).
    observation_string : :class:`str`
        Name/identifier for the observation dataset.
    weights : :class:`numpy.ndarray` or :class:`float`
        Per-vertex weights (or timeseries weights). Default is NaN (unset).
    weights_string : :class:`str`
        Name/identifier for the weights dataset.
    datatime : :class:`float`
        Time in years from start associated with the data (used to pick the
        comparison time). Default is 0.0.
    """

    # Define class string for registry and marshalling
    @classmethod
    def issm_enum_string(cls) -> str:
        return "Cfsurfacesquare" 

    # Initialise with default parameters    
    def __init__(self, other=None):
        self.name = ""
        self.definitionstring = ""
        self.surfaceid = 1
        self.model_string = ""
        self.observation = np.nan
        self.observation_string = ""
        self.weights = np.nan
        self.weights_string = ""
        self.datatime = 0.0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = "   surfacesquare:\n"
        s += '{}\n'.format(class_utils._field_display(self, 'name', 'identifier for this surfacesquare response'))
        s += '{}\n'.format(class_utils._field_display(self, 'definitionstring', 'unique output definition string (e.g. Outputdefinition1)'))
        s += '{}\n'.format(class_utils._field_display(self, 'surfaceid', '1: surface, 2: base'))
        s += '{}\n'.format(class_utils._field_display(self, 'model_string', 'string for field that is modeled'))
        s += '{}\n'.format(class_utils._field_display(self, 'observation', 'observed field compared against the model'))
        s += '{}\n'.format(class_utils._field_display(self, 'observation_string', 'string identifying observed field'))
        s += '{}\n'.format(class_utils._field_display(self, 'weights', 'weights (at vertices) applied to the misfit'))
        s += '{}\n'.format(class_utils._field_display(self, 'weights_string', 'string identifying weights'))
        s += '{}\n'.format(class_utils._field_display(self, 'datatime', 'time (years from start) for data-model comparison'))
        return s

    # Define class string
    def __str__(self):
        return "ISSM - surfacesquare Class"

    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [cf.surfacesquare] fields to 3D
        """

        # Treat "unset" as scalar NaN or arrays that are all-NaN
        if np.size(self.weights) > 1:
            if not np.all(np.isnan(self.weights)):
                self.weights = mesh._project_3d(md, vector=self.weights, type="node")
        else:
            if not np.isnan(self.weights):
                self.weights = mesh._project_3d(md, vector=self.weights, type="node")

        if np.size(self.observation) > 1:
            if not np.all(np.isnan(self.observation)):
                self.observation = mesh._project_3d(md, vector=self.observation, type="node")
        else:
            if not np.isnan(self.observation):
                self.observation = mesh._project_3d(md, vector=self.observation, type="node")

        return self

    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [cf.surfacesquare] parameters.

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

        # Basic pythonic type checks
        if not isinstance(self.name, str):
            raise TypeError("surfacesquare: 'name' must be a string")

        # MATLAB allowed Outputdefinition1..2000
        outputdef_allowed = [f"Outputdefinition{i}" for i in range(1, 2001)]

        class_utils.check_field(md, fieldname="surfacesquare.surfaceid", field=self.surfaceid, values=[1, 2])
        class_utils.check_field(md, fieldname="surfacesquare.definitionstring", field=self.definitionstring, values=outputdef_allowed)

        # observation/weights: timeseries, allow NaN/Inf (MATLAB: 'timeseries',1,'NaN',1,'Inf',1)
        class_utils.check_field(
            md,
            fieldname="surfacesquare.observation",
            field=self.observation,
            timeseries=True,
            allow_nan=True,
            allow_inf=True,
        )
        class_utils.check_field(
            md,
            fieldname="surfacesquare.weights",
            field=self.weights,
            timeseries=True,
            allow_nan=True,
            allow_inf=True,
        )

        # datatime must be within simulation window
        class_utils.check_field(
            md,
            fieldname="surfacesquare.datatime",
            field=self.datatime,
            le=md.timestepping.final_time,
        )

        return md

    # Marshall method for saving the cf.surfacesquare parameters
    def marshall_class(self, fid, prefix, md=None):
        """
        Marshall [cf.surfacesquare] parameters to a binary file.

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

        nv = md.mesh.numberofvertices
        yts = md.constants.yts

        execute._write_model_field(fid, prefix, data=self.name, name="md.cfsurfacesquare.name", format="String")
        execute._write_model_field(fid, prefix, data=self.definitionstring, name="md.cfsurfacesquare.definitionstring", format="String")
        execute._write_model_field(fid, prefix, data=self.surfaceid, name="md.cfsurfacesquare.surfaceid", format="Integer")
        execute._write_model_field(fid, prefix, data=self.model_string, name="md.cfsurfacesquare.model_string", format="String")

        execute._write_model_field(
            fid,
            prefix,
            data=self.observation,
            name="md.cfsurfacesquare.observation",
            format="DoubleMat",
            mattype=1,
            timeserieslength=nv + 1,
            yts=yts,
        )
        execute._write_model_field(fid, prefix, data=self.observation_string, name="md.cfsurfacesquare.observation_string", format="String")

        execute._write_model_field(
            fid,
            prefix,
            data=self.weights,
            name="md.cfsurfacesquare.weights",
            format="DoubleMat",
            mattype=1,
            timeserieslength=nv + 1,
            yts=yts,
        )
        execute._write_model_field(fid, prefix, data=self.weights_string, name="md.cfsurfacesquare.weights_string", format="String")

        # MATLAB: round(datatime * yts) stored as Double
        execute._write_model_field(
            fid,
            prefix,
            data=float(np.round(self.datatime * yts)),
            name="md.cfsurfacesquare.datatime",
            format="Double",
        )

## ------------------------------------------------------
## cf.surfacesquaretransient
## ------------------------------------------------------
@class_registry.register_class
class surfacesquaretransient(class_registry.manage_state):
    """
    Transient surface-square misfit response definition class for ISSM.

    This response compares a modeled field (given by ``model_string``) against a set
    of *time-series observations* using a squared-misfit form, with optional
    per-vertex weights. In ISSM convention, transient time series are typically
    stored in a matrix/vector with a time row/column included, and are marshalled
    with ``timeserieslength = numberofvertices + 1`` and ``yts`` conversion.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    name : :class:`str`
        Identifier for this response definition.
    definitionstring : :class:`str`
        Unique output definition identifier, typically ``"OutputdefinitionN"``.
        MATLAB allows ``Outputdefinition1``..``Outputdefinition2000``.
    model_string : :class:`str`
        Name of the modeled field to compare (string key used by ISSM).
        Example: ``"Surface"``, ``"Vx"``, ``"Vy"``.
    observations : :class:`numpy.ndarray` or :class:`float`
        Observed time series to compare against. By MATLAB convention this is a
        time-series array marshalled with ``timeserieslength = nv + 1``.
        Default is NaN (unset).
    weights : :class:`numpy.ndarray` or :class:`float`
        Weights applied to the misfit, typically a time-series array with the same
        convention as ``observations``. Default is NaN (unset).

    Examples
    --------
    .. code-block:: python

        >>> cf = surfacesquaretransient()
        >>> cf.name = "SurfaceAltimetry"
        >>> cf.definitionstring = "Outputdefinition1"
        >>> cf.model_string = "Surface"
        >>> cf.observations = np.vstack([md.geometry.surface, [0.0]])  # example
        >>> cf.weights = np.ones((md.mesh.numberofvertices + 1, 1))
        >>> md.outputdefinition.definitions.append(cf)
    """

    # Define class string for registry and marshalling
    @classmethod
    def issm_enum_string(cls) -> str:
        return "Cfsurfacesquaretransient"
    
    # Initialise with default parameters
    def __init__(self, other=None):
        self.name = ""
        self.definitionstring = ""
        self.model_string = ""
        self.observations = np.nan
        self.weights = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = "   surfacesquaretransient:\n"
        s += '{}\n'.format(class_utils._field_display(self, 'name', 'identifier for this surfacesquaretransient response'))
        s += '{}\n'.format(class_utils._field_display(self, 'definitionstring', 'unique output definition string (e.g. Outputdefinition1)'))
        s += '{}\n'.format(class_utils._field_display(self, 'model_string', 'string for field that is modeled'))
        s += '{}\n'.format(class_utils._field_display(self, 'observations', 'observed field time series compared against the model'))
        s += '{}\n'.format(class_utils._field_display(self, 'weights', 'weights applied to the transient square misfit'))
        return s

    # Define class string
    def __str__(self):
        return "ISSM - surfacesquaretransient Class"

    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [cf.surfacesquaretransient] fields to 3D
        """
        
        # weights
        if np.size(self.weights) > 1:
            if not np.all(np.isnan(self.weights)):
                self.weights = mesh._project_3d(md, vector=self.weights, type="node")
        else:
            if not np.isnan(self.weights):
                self.weights = mesh._project_3d(md, vector=self.weights, type="node")

        # observations
        if np.size(self.observations) > 1:
            if not np.all(np.isnan(self.observations)):
                self.observations = mesh._project_3d(md, vector=self.observations, type="node")
        else:
            if not np.isnan(self.observations):
                self.observations = mesh._project_3d(md, vector=self.observations, type="node")

        return self

    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [cf.surfacesquaretransient] parameters.

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

        if not isinstance(self.name, str):
            raise TypeError("surfacesquaretransient: 'name' must be a string")

        # MATLAB allowed Outputdefinition1..2000
        outputdef_allowed = [f"Outputdefinition{i}" for i in range(1, 2001)]

        class_utils.check_field(
            md,
            fieldname="surfacesquaretransient.definitionstring",
            field=self.definitionstring,
            values=outputdef_allowed,
        )

        # MATLAB:
        # observations size: [nv+1 NaN], NaN allowed, Inf allowed, timeseries=1
        # weights size:      [nv+1 NaN], NaN allowed, Inf allowed
        nv = md.mesh.numberofvertices

        class_utils.check_field(
            md,
            fieldname="surfacesquaretransient.observations",
            field=self.observations,
            size=(nv + 1, np.nan),
            timeseries=True,
            allow_nan=True,
            allow_inf=True,
        )

        class_utils.check_field(
            md,
            fieldname="surfacesquaretransient.weights",
            field=self.weights,
            size=(nv + 1, np.nan),
            allow_nan=True,
            allow_inf=True,
        )

        return md

    def marshall_class(self, fid, prefix, md=None):
        """
        Marshall [cf.surfacesquaretransient] parameters to a binary file.

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

        nv = md.mesh.numberofvertices
        yts = md.constants.yts

        execute._write_model_field(fid, prefix, data=self.name, name="md.cfsurfacesquaretransient.name", format="String")
        execute._write_model_field(
            fid,
            prefix,
            data=self.definitionstring,
            name="md.cfsurfacesquaretransient.definitionstring",
            format="String",
        )
        execute._write_model_field(
            fid,
            prefix,
            data=self.model_string,
            name="md.cfsurfacesquaretransient.model_string",
            format="String",
        )

        execute._write_model_field(
            fid,
            prefix,
            data=self.observations,
            name="md.cfsurfacesquaretransient.observations",
            format="DoubleMat",
            mattype=1,
            timeserieslength=nv + 1,
            yts=yts,
        )
        execute._write_model_field(
            fid,
            prefix,
            data=self.weights,
            name="md.cfsurfacesquaretransient.weights",
            format="DoubleMat",
            mattype=1,
            timeserieslength=nv + 1,
            yts=yts,
        )

## ------------------------------------------------------
## cf.surfacelogvel
## ------------------------------------------------------
@class_registry.register_class
class surfacelogvel(class_registry.manage_state):
    """
    Surface log-velocity misfit response definition class for ISSM.

    This response is typically used in transient inversions/calibration workflows
    where the misfit compares the *logarithm of surface speed* implied by modeled
    surface velocities against observed surface velocities, with optional
    per-vertex weights and an associated observation time.

    The class stores:
      - observed Vx and Vy (often in ISSM time-series format)
      - optional strings identifying those observations
      - per-vertex weights (also often a time series)
      - a ``datatime`` (years from start) telling ISSM which time the data
        corresponds to

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    name : :class:`str`
        Identifier for this response definition.
    definitionstring : :class:`str`
        Unique output definition identifier, typically ``"OutputdefinitionN"``.
        MATLAB allows ``Outputdefinition1``..``Outputdefinition2000``.
    vxobs : :class:`numpy.ndarray` or :class:`float`
        Observed x-velocity component. MATLAB marshals this as a DoubleMat time
        series with ``timeserieslength = nv + 1`` and applies a scale of
        ``1 / yts``.
        Default is NaN (unset).
    vxobs_string : :class:`str`
        Identifier/name for the observed x-velocity dataset.
    vyobs : :class:`numpy.ndarray` or :class:`float`
        Observed y-velocity component. Same marshalling conventions as ``vxobs``.
        Default is NaN (unset).
    vyobs_string : :class:`str`
        Identifier/name for the observed y-velocity dataset.
    weights : :class:`numpy.ndarray` or :class:`float`
        Weight coefficients (typically per-vertex, often time-series formatted).
        Default is NaN (unset).
    weights_string : :class:`str`
        Identifier/name for the weights dataset.
    datatime : :class:`float`
        Time in years from start associated with the data. MATLAB writes this as
        ``round(datatime * yts)`` to the binary file.

    Notes
    -----
    - The MATLAB implementation only explicitly projects ``vxobs`` (not ``vyobs``)
      during extrusion. For strict parity, we replicate that behavior here.
      If you prefer symmetry, you can project both ``vxobs`` and ``vyobs``.
    - ``vxobs`` and ``vyobs`` are written with a ``scale = 1/yts`` in MATLAB,
      matching ISSM’s internal unit conventions (m/s vs m/yr).
    """

    # Define class string for registry and marshalling
    @classmethod
    def issm_enum_string(cls) -> str:
        return "Cfsurfacelogvel"
    
    # Initialise with default parameters
    def __init__(self, other=None):
        self.name = ""
        self.definitionstring = ""
        self.vxobs = np.nan
        self.vxobs_string = ""
        self.vyobs = np.nan
        self.vyobs_string = ""
        self.weights = np.nan
        self.weights_string = ""
        self.datatime = 0.0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = "   surfacelogvel:\n"
        s += '{}\n'.format(class_utils._field_display(self, 'name', 'identifier for this surfacelogvel response'))
        s += '{}\n'.format(class_utils._field_display(self, 'definitionstring', 'unique output definition string (e.g. Outputdefinition1)'))
        s += '{}\n'.format(class_utils._field_display(self, 'vxobs', 'observed Vx used for misfit'))
        s += '{}\n'.format(class_utils._field_display(self, 'vxobs_string', 'identifier for observed Vx'))
        s += '{}\n'.format(class_utils._field_display(self, 'vyobs', 'observed Vy used for misfit'))
        s += '{}\n'.format(class_utils._field_display(self, 'vyobs_string', 'identifier for observed Vy'))
        s += '{}\n'.format(class_utils._field_display(self, 'weights', 'weights applied to the misfit'))
        s += '{}\n'.format(class_utils._field_display(self, 'weights_string', 'identifier for weights'))
        s += '{}\n'.format(class_utils._field_display(self, 'datatime', 'time (years from start) for data-model comparison'))
        return s

    # Define class string
    def __str__(self):
        return "ISSM - surfacelogvel Class"

    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude [cf.surfacelogvel] fields to 3D
        """

        # weights
        if np.size(self.weights) > 1:
            if not np.all(np.isnan(self.weights)):
                self.weights = mesh._project_3d(md, vector=self.weights, type="node")
        else:
            if not np.isnan(self.weights):
                self.weights = mesh._project_3d(md, vector=self.weights, type="node")

        # vxobs (MATLAB projects only vxobs)
        if np.size(self.vxobs) > 1:
            if not np.all(np.isnan(self.vxobs)):
                self.vxobs = mesh._project_3d(md, vector=self.vxobs, type="node")
        else:
            if not np.isnan(self.vxobs):
                self.vxobs = mesh._project_3d(md, vector=self.vxobs, type="node")

        return self

    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [cf.surfacelogvel] parameters.

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

        if not isinstance(self.name, str):
            raise TypeError("surfacelogvel: 'name' must be a string")

        # MATLAB allowed Outputdefinition1..2000
        outputdef_allowed = [f"Outputdefinition{i}" for i in range(1, 2001)]

        class_utils.check_field(
            md,
            fieldname="surfacelogvel.definitionstring",
            field=self.definitionstring,
            values=outputdef_allowed,
        )

        # MATLAB checks only vxobs, weights, datatime (vyobs is not checked there).
        # We keep parity, but you can add vyobs similarly if desired.
        class_utils.check_field(
            md,
            fieldname="surfacelogvel.vxobs",
            field=self.vxobs,
            timeseries=True,
            allow_nan=True,
            allow_inf=True,
        )
        class_utils.check_field(
            md,
            fieldname="surfacelogvel.weights",
            field=self.weights,
            timeseries=True,
            allow_nan=True,
            allow_inf=True,
        )
        class_utils.check_field(
            md,
            fieldname="surfacelogvel.datatime",
            field=self.datatime,
            le=md.timestepping.final_time,
        )

        return md

    # Marshall method for saving the cf.surfacelogvel parameters
    def marshall_class(self, fid, prefix, md=None):
        """
        Marshall [cf.surfacelogvel] parameters to a binary file.

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

        nv = md.mesh.numberofvertices
        yts = md.constants.yts

        execute._write_model_field(fid, prefix, data=self.name, name="md.cfsurfacelogvel.name", format="String")
        execute._write_model_field(
            fid,
            prefix,
            data=self.definitionstring,
            name="md.cfsurfacelogvel.definitionstring",
            format="String",
        )

        execute._write_model_field(
            fid,
            prefix,
            data=self.vxobs,
            name="md.cfsurfacelogvel.vxobs",
            format="DoubleMat",
            mattype=1,
            timeserieslength=nv + 1,
            yts=yts,
            scale=1.0 / yts,
        )
        execute._write_model_field(
            fid,
            prefix,
            data=self.vxobs_string,
            name="md.cfsurfacelogvel.vxobs_string",
            format="String",
        )

        execute._write_model_field(
            fid,
            prefix,
            data=self.vyobs,
            name="md.cfsurfacelogvel.vyobs",
            format="DoubleMat",
            mattype=1,
            timeserieslength=nv + 1,
            yts=yts,
            scale=1.0 / yts,
        )
        execute._write_model_field(
            fid,
            prefix,
            data=self.vyobs_string,
            name="md.cfsurfacelogvel.vyobs_string",
            format="String",
        )

        execute._write_model_field(
            fid,
            prefix,
            data=self.weights,
            name="md.cfsurfacelogvel.weights",
            format="DoubleMat",
            mattype=1,
            timeserieslength=nv + 1,
            yts=yts,
        )
        execute._write_model_field(
            fid,
            prefix,
            data=self.weights_string,
            name="md.cfsurfacelogvel.weights_string",
            format="String",
        )

        execute._write_model_field(
            fid,
            prefix,
            data=float(np.round(self.datatime * yts)),
            name="md.cfsurfacelogvel.datatime",
            format="Double",
        )
