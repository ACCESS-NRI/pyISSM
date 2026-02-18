import numpy as np

from pyissm.model.classes import class_utils, class_registry
from pyissm.model import execute, mesh


@class_registry.register_class
class cflevelsetmisfit(class_registry.manage_state):
    """
    Level-set misfit response definition.

    Python/pyISSM equivalent of ISSM MATLAB's ``cflevelsetmisfit`` class.

    This response is commonly used to measure misfit between a modeled level-set
    field (e.g. calving-front / ice-mask level set) and an observed level-set
    field. It supports optional per-vertex weights and an associated data time
    (years from start) indicating when the observation applies.

    Parameters
    ----------
    other : object, optional
        Any object with matching attribute names. If provided, its values overwrite
        the defaults (pyISSM convention).

    Attributes
    ----------
    name : str
        Identifier for this response definition.
    definitionstring : str
        Unique output definition identifier, typically ``"OutputdefinitionN"``.
        MATLAB allows ``Outputdefinition1``..``Outputdefinition2000``.
    model_string : str
        String key for the modeled field to be compared.
        Example: ``"MaskIceLevelset"``.
    observation : array_like or float
        Observed field to compare against (often in ISSM time-series format).
        Written as a DoubleMat time series with ``timeserieslength = nv + 1``.
        Default is NaN (unset).
    observation_string : str
        Identifier/name for the observation dataset.
    weights : array_like or float
        Weight coefficients (per-vertex, often time-series formatted).
        Written as a DoubleMat time series with ``timeserieslength = nv + 1``.
        Default is NaN (unset).
    weights_string : str
        Identifier/name for the weights dataset.
    datatime : float
        Time in years from start associated with the data. MATLAB writes this as
        ``round(datatime * yts)`` to the binary file.

    Methods
    -------
    extrude(md)
        Project node-based fields to 3D (no-op if unset).
    check_consistency(md, solution=None, analyses=None)
        Validate field types and ranges.
    marshall_class(fid, prefix, md=None)
        Write this definition to ISSM binary format via ``execute.WriteData``.

    Examples
    --------
    >>> cf = cflevelsetmisfit()
    >>> cf.name = "CalvingFrontPosition"
    >>> cf.definitionstring = "Outputdefinition1"
    >>> cf.model_string = "MaskIceLevelset"
    >>> cf.observation_string = "LevelsetObservations"
    >>> cf.observation = md.mask.ice_levelset
    >>> cf.weights = np.ones((md.mesh.numberofvertices, 1))
    >>> cf.weights_string = "WeightsLevelsetObservations"
    >>> cf.datatime = time
    >>> md.outputdefinition.definitions.append(cf)
    """

    def __init__(self, other=None):
        # Defaults (MATLAB parity)
        self.name = ""
        self.definitionstring = ""
        self.model_string = ""
        self.observation = np.nan
        self.observation_string = ""
        self.weights = np.nan
        self.weights_string = ""
        self.datatime = 0.0

        super().__init__(other)

    def __repr__(self):
        s = "   cflevelsetmisfit:\n"
        s += f"{class_utils.fielddisplay(self, 'name', 'identifier for this cflevelsetmisfit response')}\n"
        s += f"{class_utils.fielddisplay(self, 'definitionstring', 'unique output definition string (e.g. Outputdefinition1)')}\n"
        s += f"{class_utils.fielddisplay(self, 'model_string', 'string for field that is modeled')}\n"
        s += f"{class_utils.fielddisplay(self, 'observation', 'observed field compared against the model')}\n"
        s += f"{class_utils.fielddisplay(self, 'observation_string', 'identifier for observed field')}\n"
        s += f"{class_utils.fielddisplay(self, 'weights', 'weights applied to the misfit')}\n"
        s += f"{class_utils.fielddisplay(self, 'weights_string', 'identifier for weights')}\n"
        s += f"{class_utils.fielddisplay(self, 'datatime', 'time (years from start) for data-model comparison')}\n"
        return s

    def __str__(self):
        return "ISSM - cflevelsetmisfit Class"

    def extrude(self, md):
        """
        Extrude node-based fields to 3D.

        Projects ``weights`` and ``observation`` with ``mesh.project_3d`` if they
        are set (i.e., not entirely NaN). Mirrors MATLAB behavior.
        """
        # weights
        if np.size(self.weights) > 1:
            if not np.all(np.isnan(self.weights)):
                self.weights = mesh.project_3d(md, vector=self.weights, type="node")
        else:
            if not np.isnan(self.weights):
                self.weights = mesh.project_3d(md, vector=self.weights, type="node")

        # observation
        if np.size(self.observation) > 1:
            if not np.all(np.isnan(self.observation)):
                self.observation = mesh.project_3d(md, vector=self.observation, type="node")
        else:
            if not np.isnan(self.observation):
                self.observation = mesh.project_3d(md, vector=self.observation, type="node")

        return self

    def check_consistency(self, md, solution=None, analyses=None):
        """
        Check field validity for this response definition.
        """
        if not isinstance(self.name, str):
            raise TypeError("cflevelsetmisfit: 'name' must be a string")

        # MATLAB allowed Outputdefinition1..2000
        outputdef_allowed = [f"Outputdefinition{i}" for i in range(1, 2001)]

        class_utils.check_field(
            md,
            fieldname="cflevelsetmisfit.definitionstring",
            field=self.definitionstring,
            values=outputdef_allowed,
        )

        class_utils.check_field(
            md,
            fieldname="cflevelsetmisfit.observation",
            field=self.observation,
            timeseries=True,
            allow_nan=True,
            allow_inf=True,
        )
        class_utils.check_field(
            md,
            fieldname="cflevelsetmisfit.weights",
            field=self.weights,
            timeseries=True,
            allow_nan=True,
            allow_inf=True,
        )
        class_utils.check_field(
            md,
            fieldname="cflevelsetmisfit.datatime",
            field=self.datatime,
            le=md.timestepping.final_time,
        )

        return md

    def marshall_class(self, fid, prefix, md=None):
        """
        Marshall this response definition to ISSM binary input.

        Mirrors MATLAB marshalling:
          - observation/weights written as DoubleMat time series with:
              mattype=1, timeserieslength=nv+1, yts conversion
          - datatime written as ``round(datatime * yts)`` (seconds) as Double
        """
        nv = md.mesh.numberofvertices
        yts = md.constants.yts

        execute.WriteData(fid, prefix, data=self.name, name="md.cflevelsetmisfit.name", format="String")
        execute.WriteData(
            fid,
            prefix,
            data=self.definitionstring,
            name="md.cflevelsetmisfit.definitionstring",
            format="String",
        )
        execute.WriteData(
            fid,
            prefix,
            data=self.model_string,
            name="md.cflevelsetmisfit.model_string",
            format="String",
        )

        execute.WriteData(
            fid,
            prefix,
            data=self.observation,
            name="md.cflevelsetmisfit.observation",
            format="DoubleMat",
            mattype=1,
            timeserieslength=nv + 1,
            yts=yts,
        )
        execute.WriteData(
            fid,
            prefix,
            data=self.observation_string,
            name="md.cflevelsetmisfit.observation_string",
            format="String",
        )

        execute.WriteData(
            fid,
            prefix,
            data=self.weights,
            name="md.cflevelsetmisfit.weights",
            format="DoubleMat",
            mattype=1,
            timeserieslength=nv + 1,
            yts=yts,
        )
        execute.WriteData(
            fid,
            prefix,
            data=self.weights_string,
            name="md.cflevelsetmisfit.weights_string",
            format="String",
        )

        execute.WriteData(
            fid,
            prefix,
            data=float(np.round(self.datatime * yts)),
            name="md.cflevelsetmisfit.datatime",
            format="Double",
        )
