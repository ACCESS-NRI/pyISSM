import numpy as np
import warnings

from pyissm.model.classes import class_utils, class_registry
from pyissm.model import execute, mesh


@class_registry.register_class
class cflevelsetmisfit(class_registry.manage_state):
    """
    MISFIT: Levelset (calving front position) misfit definition.

    MATLAB usage:
      cflevelsetmisfit('name','CalvingFrontPosition',
                       'definitionstring','Outputdefinition1',
                       'model_string','MaskIceLevelset',
                       'observation_string','LevelsetObservations',
                       'observation',md.mask.ice_levelset,
                       'weights',ones(md.mesh.numberofvertices,1),
                       'weights_string','WeightsLevelsetObservations',
                       'datatime',time);
    """

    def __init__(self, other=None):
        # fields
        self.name = ''
        self.definitionstring = ''     # 'Outputdefinition#'
        self.model_string = ''         # modeled field name (e.g. 'MaskIceLevelset')
        self.observation = np.nan      # observed field (vertex, possibly timeseries)
        self.observation_string = ''   # identifier for observation
        self.weights = np.nan          # weights (vertex, possibly timeseries)
        self.weights_string = ''       # identifier for weights
        self.datatime = 0.0            # time in years from start

        # inherit from provided class/struct-like object
        super().__init__(other)

    def __repr__(self):
        s = "   cflevelsetmisfit:\n"
        s += f"{class_utils.fielddisplay(self, 'name', 'identifier for this cflevelsetmisfit response')}\n"
        s += f"{class_utils.fielddisplay(self, 'definitionstring', 'unique identifier string Outputdefinition#')}\n"
        s += f"{class_utils.fielddisplay(self, 'model_string', 'string for field that is modeled')}\n"
        s += f"{class_utils.fielddisplay(self, 'observation', 'observed field that we compare the model against')}\n"
        s += f"{class_utils.fielddisplay(self, 'observation_string', 'observation string')}\n"
        s += f"{class_utils.fielddisplay(self, 'weights', 'weights (at vertices) to apply to the cflevelsetmisfit')}\n"
        s += f"{class_utils.fielddisplay(self, 'weights_string', 'string for weights for identification purposes')}\n"
        s += f"{class_utils.fielddisplay(self, 'datatime', 'time to compare data to model for misfit')}\n"
        return s

    def __str__(self):
        return "ISSM - cflevelsetmisfit Class"

    def extrude(self, md):
        """
        Extrude observation/weights to 3D (node-based), matching MATLAB.
        """
        if not np.isscalar(self.weights) and not (isinstance(self.weights, float) and np.isnan(self.weights)):
            self.weights = mesh.project_3d(md, vector=self.weights, type='node')
        if not np.isscalar(self.observation) and not (isinstance(self.observation, float) and np.isnan(self.observation)):
            self.observation = mesh.project_3d(md, vector=self.observation, type='node')
        return self

    def check_consistency(self, md, solution, analyses):
        if not isinstance(self.name, str):
            raise TypeError("cflevelsetmisfit: 'name' field should be a string!")

        # Definitionstring must look like Outputdefinition#
        # (MATLAB allows up to 2000; keeping same spirit)
        allowed = [f"Outputdefinition{i}" for i in range(1, 2001)]
        class_utils.check_field(
            md,
            fieldname="self.definitionstring",
            field=self.definitionstring,
            values=allowed,
        )

        # observation/weights are vertex timeseries (NaN/Inf allowed per MATLAB)
        class_utils.check_field(
            md,
            fieldname="self.observation",
            field=self.observation,
            timeseries=True,
            allow_nan=True,
            allow_inf=True,
        )
        class_utils.check_field(
            md,
            fieldname="self.weights",
            field=self.weights,
            timeseries=True,
            allow_nan=True,
            allow_inf=True,
        )

        # datatime <= final_time (MATLAB)
        class_utils.check_field(
            md,
            fieldname="self.datatime",
            field=self.datatime,
            le=md.timestepping.final_time,
        )
        return md

    def marshall_class(self, fid, prefix, md=None, write_type=True):
        """
        Marshall to binary, matching MATLAB marshall().
        """
        # Strings
        execute.WriteData(fid, prefix, data=self.name,
                          name="md.cflevelsetmisfit.name", format="String")
        execute.WriteData(fid, prefix, data=self.definitionstring,
                          name="md.cflevelsetmisfit.definitionstring", format="String")
        execute.WriteData(fid, prefix, data=self.model_string,
                          name="md.cflevelsetmisfit.model_string", format="String")

        # observation (vertex timeseries)
        execute.WriteData(
            fid, prefix,
            data=self.observation,
            name="md.cflevelsetmisfit.observation",
            format="DoubleMat",
            mattype=1,
            timeserieslength=md.mesh.numberofvertices + 1,
            yts=md.constants.yts,
        )

        # observation string
        execute.WriteData(fid, prefix, data=self.observation_string,
                          name="md.cflevelsetmisfit.observation_string", format="String")

        # weights (vertex timeseries)
        execute.WriteData(
            fid, prefix,
            data=self.weights,
            name="md.cflevelsetmisfit.weights",
            format="DoubleMat",
            mattype=1,
            timeserieslength=md.mesh.numberofvertices + 1,
            yts=md.constants.yts,
        )

        # weights string
        execute.WriteData(fid, prefix, data=self.weights_string,
                          name="md.cflevelsetmisfit.weights_string", format="String")

        # datatime stored in seconds (rounded), matching MATLAB:
        # round(self.datatime*md.constants.yts)
        execute.WriteData(
            fid, prefix,
            data=float(np.round(self.datatime * md.constants.yts)),
            name="md.cflevelsetmisfit.datatime",
            format="Double",
        )
