from __future__ import annotations
from dataclasses import dataclass
from typing import Any
import numpy as np

from pyissm.model.mesh import project_3d
from pyissm.model.classes.param_utils import check_field
from pyissm.model.io import WriteData


def _is_unset(x: Any) -> bool:
    """Match MATLAB-ish semantics: NaN means 'unset' (scalar or array)."""
    if x is None:
        return True
    if isinstance(x, (float, np.floating)):
        return bool(np.isnan(x))
    try:
        arr = np.asarray(x)
        if arr.size == 0:
            return True
        if np.issubdtype(arr.dtype, np.number):
            return bool(np.all(np.isnan(arr)))
    except Exception:
        pass
    return False


@dataclass
class CFSurfaceLogVel:
    # cfsurfacelogvel
    name: str = ""
    definitionstring: str = ""  # 'Outputdefinition[1-2000]'
    vxobs: Any = np.nan
    vxobs_string: str = ""
    vyobs: Any = np.nan
    vyobs_string: str = ""
    weights: Any = np.nan
    weights_string: str = ""
    datatime: float = 0.0  # years from start

    def setdefaultparameters(self) -> "CFSurfaceLogVel":
        self.datatime = 0.0
        return self

    def extrude(self, md) -> "CFSurfaceLogVel":
        if not _is_unset(self.weights):
            self.weights = project_3d(md, "vector", self.weights, "type", "node")
        if not _is_unset(self.vxobs):
            self.vxobs = project_3d(md, "vector", self.vxobs, "type", "node")
        # MATLAB version did NOT project vyobs; keep identical behavior.
        return self

    def checkconsistency(self, md, solution=None, analyses=None):
        if not isinstance(self.name, str):
            raise TypeError("cfsurfacelogvel: 'name' field should be a string!")

        allowed = [f"Outputdefinition{i}" for i in range(1, 2001)]
        md = check_field(
            md,
            fieldname="self.definitionstring",
            field=self.definitionstring,
            values=allowed,
        )

        md = check_field(
            md,
            fieldname="self.vxobs",
            field=self.vxobs,
            timeseries=1,
            NaN=1,
            Inf=1,
        )
        md = check_field(
            md,
            fieldname="self.weights",
            field=self.weights,
            timeseries=1,
            NaN=1,
            Inf=1,
        )

        md = check_field(
            md,
            fieldname="self.datatime",
            field=self.datatime,
            **{"<=": md.timestepping.final_time},  # mirrors MATLAB checkfield(...,'<=',...)
        )

        return md

    def marshall(self, prefix: str, md, fid):
        yts = md.constants.yts
        nverts = md.mesh.numberofvertices

        WriteData(fid, prefix, data=self.name,
                  name="md.cfsurfacelogvel.name", format="String")
        WriteData(fid, prefix, data=self.definitionstring,
                  name="md.cfsurfacelogvel.definitionstring", format="String")

        WriteData(
            fid, prefix,
            data=self.vxobs,
            name="md.cfsurfacelogvel.vxobs",
            format="DoubleMat",
            mattype=1,
            timeserieslength=nverts + 1,
            yts=yts,
            scale=1.0 / yts,
        )
        WriteData(fid, prefix, data=self.vxobs_string,
                  name="md.cfsurfacelogvel.vxobs_string", format="String")

        WriteData(
            fid, prefix,
            data=self.vyobs,
            name="md.cfsurfacelogvel.vyobs",
            format="DoubleMat",
            mattype=1,
            timeserieslength=nverts + 1,
            yts=yts,
            scale=1.0 / yts,
        )
        WriteData(fid, prefix, data=self.vyobs_string,
                  name="md.cfsurfacelogvel.vyobs_string", format="String")

        WriteData(
            fid, prefix,
            data=self.weights,
            name="md.cfsurfacelogvel.weights",
            format="DoubleMat",
            mattype=1,
            timeserieslength=nverts + 1,
            yts=yts,
        )
        WriteData(fid, prefix, data=self.weights_string,
                  name="md.cfsurfacelogvel.weights_string", format="String")

        WriteData(
            fid, prefix,
            data=float(np.round(self.datatime * yts)),
            name="md.cfsurfacelogvel.datatime",
            format="Double",
        )
