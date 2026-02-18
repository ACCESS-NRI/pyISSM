import numpy as np
from pyissm.model.classes import class_utils, class_registry
from pyissm.model import execute, mesh


@class_registry.register_class
class cflevelsetmisfit(class_registry.manage_state):
    """
    Levelset misfit (calving front position).

    Example usage:

        cf = cflevelsetmisfit(
            name="CalvingFrontPosition",
            definitionstring="Outputdefinition1",
            model_string="MaskIceLevelset",
            observation_string="LevelsetObservations",
            observation=md.mask.ice_levelset,
            weights=np.ones(md.mesh.numberofvertices),
            weights_string="WeightsLevelsetObservations",
            datatime=time,
        )
    """

    def __init__(self, other=None, **kwargs):
        # --- defaults ---
        self.name = ""
        self.definitionstring = ""
        self.model_string = ""
        self.observation = np.nan
        self.observation_string = ""
        self.weights = np.nan
        self.weights_string = ""
        self.datatime = 0.0

        # inherit from another instance if provided
        super().__init__(other)

        # apply keyword overrides
        for key, value in kwargs.items():
            if not hasattr(self, key):
                raise AttributeError(f"cflevelsetmisfit: unknown field '{key}'")
            setattr(self, key, value)

    # -------------------------------------------------
    # Extrude (node-based)
    # -------------------------------------------------
    def extrude(self, md):
        if not np.isscalar(self.weights):
            self.weights = mesh.project_3d(md, vector=self.weights, type="node")
        if not np.isscalar(self.observation):
            self.observation = mesh.project_3d(md, vector=self.observation, type="node")
        return self

    # -------------------------------------------------
    # Consistency
    # -------------------------------------------------
    def check_consistency(self, md, solution, analyses):

        if not isinstance(self.name, str):
            raise TypeError("cflevelsetmisfit: 'name' must be a string")

        allowed = [f"Outputdefinition{i}" for i in range(1, 2001)]
        class_utils.check_field(
            md,
            fieldname="self.definitionstring",
            field=self.definitionstring,
            values=allowed,
        )

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

        class_utils.check_field(
            md,
            fieldname="self.datatime",
            field=self.datatime,
            le=md.timestepping.final_time,
        )

        return md

    # -------------------------------------------------
    # Marshall
    # -------------------------------------------------
    def marshall_class(self, fid, prefix, md=None, write_type=True):

        execute.WriteData(fid, prefix,
                          data=self.name,
                          name="md.cflevelsetmisfit.name",
                          format="String")

        execute.WriteData(fid, prefix,
                          data=self.definitionstring,
                          name="md.cflevelsetmisfit.definitionstring",
                          format="String")

        execute.WriteData(fid, prefix,
                          data=self.model_string,
                          name="md.cflevelsetmisfit.model_string",
                          format="String")

        execute.WriteData(
            fid, prefix,
            data=self.observation,
            name="md.cflevelsetmisfit.observation",
            format="DoubleMat",
            mattype=1,
            timeserieslength=md.mesh.numberofvertices + 1,
            yts=md.constants.yts,
        )

        execute.WriteData(fid, prefix,
                          data=self.observation_string,
                          name="md.cflevelsetmisfit.observation_string",
                          format="String")

        execute.WriteData(
            fid, prefix,
            data=self.weights,
            name="md.cflevelsetmisfit.weights",
            format="DoubleMat",
            mattype=1,
            timeserieslength=md.mesh.numberofvertices + 1,
            yts=md.constants.yts,
        )

        execute.WriteData(fid, prefix,
                          data=self.weights_string,
                          name="md.cflevelsetmisfit.weights_string",
                          format="String")

        execute.WriteData(
            fid, prefix,
            data=float(np.round(self.datatime * md.constants.yts)),
            name="md.cflevelsetmisfit.datatime",
            format="Double",
        )
