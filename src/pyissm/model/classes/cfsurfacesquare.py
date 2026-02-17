import numpy as np

from pyissm.model.classes import class_utils, class_registry
from pyissm.model import execute, mesh


@class_registry.register_class
class cfsurfacesquare(class_registry.manage_state):
    """
    Surface-square cost-function (response) definition.

    Python/pyISSM equivalent of ISSM MATLAB's ``cfsurfacesquare`` class. This
    response is typically used to compare a modeled surface (or base) field
    against observed data with per-vertex weights, optionally as a time series.

    Parameters
    ----------
    other : object, optional
        Any object with matching attributes. If provided, its values overwrite
        the defaults.

    Attributes
    ----------
    name : str
        Identifier for this response definition.
    definitionstring : str
        Unique output definition identifier, typically ``"OutputdefinitionN"``.
    surfaceid : int
        Which surface to evaluate:
        - 1: surface
        - 2: base
    model_string : str
        Name of the modeled field to compare (string key used by ISSM).
    observation : array_like or float
        Observed field to compare against. May be provided as a time series in
        ISSM “timeseries” format (length = numberofvertices + 1).
        Default is NaN (unset).
    observation_string : str
        Name/identifier for the observation dataset.
    weights : array_like or float
        Per-vertex weights (or timeseries weights). Default is NaN (unset).
    weights_string : str
        Name/identifier for the weights dataset.
    datatime : float
        Time in years from start associated with the data (used to pick the
        comparison time). Default is 0.0.

    Notes
    -----
    - The MATLAB implementation treats ``observation`` and ``weights`` as
      time series and writes them with ``timeserieslength = nv + 1`` and
      ``yts`` conversion. This class mirrors that marshalling behavior.
    """

    def __init__(self, other=None):
        # Defaults (MATLAB parity)
        self.name = ""
        self.definitionstring = ""
        self.surfaceid = 1
        self.model_string = ""
        self.observation = np.nan
        self.observation_string = ""
        self.weights = np.nan
        self.weights_string = ""
        self.datatime = 0.0

        super().__init__(other)

    def __repr__(self):
        s = "   cfsurfacesquare:\n"
        s += f"{class_utils.fielddisplay(self, 'name', 'identifier for this cfsurfacesquare response')}\n"
        s += f"{class_utils.fielddisplay(self, 'definitionstring', 'unique output definition string (e.g. Outputdefinition1)')}\n"
        s += f"{class_utils.fielddisplay(self, 'surfaceid', '1: surface, 2: base')}\n"
        s += f"{class_utils.fielddisplay(self, 'model_string', 'string for field that is modeled')}\n"
        s += f"{class_utils.fielddisplay(self, 'observation', 'observed field compared against the model')}\n"
        s += f"{class_utils.fielddisplay(self, 'observation_string', 'string identifying observed field')}\n"
        s += f"{class_utils.fielddisplay(self, 'weights', 'weights (at vertices) applied to the misfit')}\n"
        s += f"{class_utils.fielddisplay(self, 'weights_string', 'string identifying weights')}\n"
        s += f"{class_utils.fielddisplay(self, 'datatime', 'time (years from start) for data-model comparison')}\n"
        return s

    def __str__(self):
        return "ISSM - cfsurfacesquare Class"

    def extrude(self, md):
        """
        Extrude node-based fields to 3D.

        Projects ``weights`` and ``observation`` with ``mesh.project_3d`` if they
        are provided (i.e., not entirely NaN).
        """
        # Treat "unset" as scalar NaN or arrays that are all-NaN
        if np.size(self.weights) > 1:
            if not np.all(np.isnan(self.weights)):
                self.weights = mesh.project_3d(md, vector=self.weights, type="node")
        else:
            if not np.isnan(self.weights):
                self.weights = mesh.project_3d(md, vector=self.weights, type="node")

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
        # Basic pythonic type checks
        if not isinstance(self.name, str):
            raise TypeError("cfsurfacesquare: 'name' must be a string")

        # MATLAB allowed Outputdefinition1..2000
        outputdef_allowed = [f"Outputdefinition{i}" for i in range(1, 2001)]

        class_utils.check_field(md, fieldname="cfsurfacesquare.surfaceid", field=self.surfaceid, values=[1, 2])
        class_utils.check_field(md, fieldname="cfsurfacesquare.definitionstring", field=self.definitionstring, values=outputdef_allowed)

        # observation/weights: timeseries, allow NaN/Inf (MATLAB: 'timeseries',1,'NaN',1,'Inf',1)
        class_utils.check_field(
            md,
            fieldname="cfsurfacesquare.observation",
            field=self.observation,
            timeseries=True,
            allow_nan=True,
            allow_inf=True,
        )
        class_utils.check_field(
            md,
            fieldname="cfsurfacesquare.weights",
            field=self.weights,
            timeseries=True,
            allow_nan=True,
            allow_inf=True,
        )

        # datatime must be within simulation window
        class_utils.check_field(
            md,
            fieldname="cfsurfacesquare.datatime",
            field=self.datatime,
            le=md.timestepping.final_time,
        )

        return md

    def marshall_class(self, fid, prefix, md=None):
        """
        Marshall this response definition to ISSM binary input.

        Mirrors MATLAB:
        - writes strings + ints directly
        - writes observation/weights as timeseries DoubleMat with length nv+1 and yts
        - writes datatime in seconds (rounded) as Double
        """
        nv = md.mesh.numberofvertices
        yts = md.constants.yts

        execute.WriteData(fid, prefix, data=self.name, name="md.cfsurfacesquare.name", format="String")
        execute.WriteData(fid, prefix, data=self.definitionstring, name="md.cfsurfacesquare.definitionstring", format="String")
        execute.WriteData(fid, prefix, data=self.surfaceid, name="md.cfsurfacesquare.surfaceid", format="Integer")
        execute.WriteData(fid, prefix, data=self.model_string, name="md.cfsurfacesquare.model_string", format="String")

        execute.WriteData(
            fid,
            prefix,
            data=self.observation,
            name="md.cfsurfacesquare.observation",
            format="DoubleMat",
            mattype=1,
            timeserieslength=nv + 1,
            yts=yts,
        )
        execute.WriteData(fid, prefix, data=self.observation_string, name="md.cfsurfacesquare.observation_string", format="String")

        execute.WriteData(
            fid,
            prefix,
            data=self.weights,
            name="md.cfsurfacesquare.weights",
            format="DoubleMat",
            mattype=1,
            timeserieslength=nv + 1,
            yts=yts,
        )
        execute.WriteData(fid, prefix, data=self.weights_string, name="md.cfsurfacesquare.weights_string", format="String")

        # MATLAB: round(datatime * yts) stored as Double
        execute.WriteData(
            fid,
            prefix,
            data=float(np.round(self.datatime * yts)),
            name="md.cfsurfacesquare.datatime",
            format="Double",
        )
