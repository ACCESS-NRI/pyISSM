import numpy as np
import pyissm

# ---------------------------------------------------------------------
# ASSUMPTION:
# - You already ran the "Generate observations" transient solve
# - Then you reset md.materials.rheology_B[0:-1,:] = 1.8e8 (constant guess)
# - Now you want to build the cost functions + autodiff control solve
# ---------------------------------------------------------------------

# Convenience handles (you may need to adjust these import paths)
# In some pyISSM installs these live under pyissm.model.classes / pyissm.model.outputdefinition / etc.
from pyissm.model.classes.independent import independent
from pyissm.model.classes.dependent import dependent

# Cost-function outputdefinition constructors (ADJUST PATHS IF NEEDED)
# Often these are e.g. pyissm.model.outputdefinition.cfsurfacesquare, etc.
from pyissm.model.classes import (
    cfsurfacelogvel,
    cfsurfacesquare,
)

# Ensure containers exist
# (Depending on pyISSM version, these may already be initialized)
if not hasattr(md, "outputdefinition") or md.outputdefinition is None:
    md.outputdefinition = pyissm.model.classes.outputdefinition.outputdefinition()

if not hasattr(md.outputdefinition, "definitions") or md.outputdefinition.definitions is None:
    md.outputdefinition.definitions = {}

if not hasattr(md, "autodiff") or md.autodiff is None:
    md.autodiff = pyissm.model.classes.autodiff.autodiff()

# In MATLAB they do md.autodiff.dependents{count} = ...
# In pyISSM this is commonly a Python list. Ensure it exists.
if not hasattr(md.autodiff, "dependents") or md.autodiff.dependents is None:
    md.autodiff.dependents = []

# ---------------------------------------------------------------------
# Set cost function outputdefinitions + dependents (one set per time slice)
# ---------------------------------------------------------------------
count = 0  # Python 0-based; we'll append dependents

for sol in md.results.TransientSolution:
    vx_obs = sol.Vx
    vy_obs = sol.Vy
    z_obs  = sol.Surface
    time   = sol.time

    weights = np.ones((md.mesh.numberofvertices,))

    # 1) LogVel misfit
    name = f"LogVelMis{count+1}"
    defstring = f"Outputdefinition{count+1}"
    md.outputdefinition.definitions[count+1] = cfsurfacelogvel(
        name=name,
        definitionstring=defstring,
        vxobs_string="VxObs",
        vxobs=vx_obs,
        vyobs_string="VyObs",
        vyobs=vy_obs,
        weights=weights,
        weights_string="WeightsSurfaceObservation",
        datatime=time,
    )
    md.autodiff.dependents.append(
        dependent(
            name=defstring,
            type="scalar",
            fos_reverse_index=1,
        )
    )
    count += 1

    # 2) Vy square misfit
    name = f"VyMisfit{count+1}"
    defstring = f"Outputdefinition{count+1}"
    md.outputdefinition.definitions[count+1] = cfsurfacesquare(
        name=name,
        definitionstring=defstring,
        model_string="Vy",
        observation_string="VyObs",
        observation=vy_obs / md.constants.yts,
        weights=weights,
        weights_string="WeightsSurfaceObservation",
        datatime=time,
    )
    md.autodiff.dependents.append(
        dependent(name=defstring, type="scalar", fos_reverse_index=1)
    )
    count += 1

    # 3) Vx square misfit (note MATLAB uses 500*weights)
    name = f"VxMisfit{count+1}"
    defstring = f"Outputdefinition{count+1}"
    md.outputdefinition.definitions[count+1] = cfsurfacesquare(
        name=name,
        definitionstring=defstring,
        model_string="Vx",
        observation_string="VxObs",
        observation=vx_obs / md.constants.yts,
        weights=500.0 * weights,
        weights_string="WeightsSurfaceObservation",
        datatime=time,
    )
    md.autodiff.dependents.append(
        dependent(name=defstring, type="scalar", fos_reverse_index=1)
    )
    count += 1

    # 4) DEM / Surface square misfit (note MATLAB uses (1/yts)*weights)
    name = f"DEMMisfit{count+1}"
    defstring = f"Outputdefinition{count+1}"
    md.outputdefinition.definitions[count+1] = cfsurfacesquare(
        name=name,
        definitionstring=defstring,
        model_string="Surface",
        observation_string="SurfaceObservation",
        observation=z_obs,
        weights=(1.0 / md.constants.yts) * weights,
        weights_string="WeightsSurfaceObservation",
        datatime=time,
    )
    md.autodiff.dependents.append(
        dependent(name=defstring, type="scalar", fos_reverse_index=1)
    )
    count += 1


# ---------------------------------------------------------------------
# Independent (control): MaterialsRheologyBbar with bounds
# MATLAB:
# min_params = md.materials.rheology_B; min_params(1:end-1,:) = cuffey(273);
# max_params = md.materials.rheology_B; max_params(1:end-1,:) = cuffey(200);
# ---------------------------------------------------------------------

min_params = md.materials.rheology_B.copy()
max_params = md.materials.rheology_B.copy()

# NOTE:
# - "cuffey(T)" is an ISSM helper giving B(T). In pyISSM you might have:
#   pyissm.tools.materials.cuffey(temperatureK)
#   or pyissm.tools.materials.paterson() etc.
# Replace these two lines with the correct function in your build.
Bmin = pyissm.tools.materials.cuffey(273.0)  # <-- adjust if needed
Bmax = pyissm.tools.materials.cuffey(200.0)  # <-- adjust if needed

# apply to all vertices, both time columns, excluding last "time row" (0:-1)
min_params[0:-1, :] = Bmin
max_params[0:-1, :] = Bmax

md.autodiff.independents = [
    independent(
        name="MaterialsRheologyBbar",
        control_size=md.materials.rheology_B.shape[1],   # number of columns (=2)
        type="vertex",
        min_parameters=min_params,
        max_parameters=max_params,
        control_scaling_factor=1e8,
        fos_reverse_index=1,  # keep consistent with dependents if your API expects it here
    )
]

# ---------------------------------------------------------------------
# Inversion + autodiff flags
# ---------------------------------------------------------------------
md.inversion = pyissm.model.inversion.adm1qn3inversion(md.inversion)
md.inversion.iscontrol = 1
md.inversion.maxiter = 4
md.inversion.maxsteps = md.inversion.maxiter
md.inversion.dxmin = 1e-5

md.autodiff.isautodiff = True
md.autodiff.driver = "fos_reverse"

md.verbose = pyissm.model.verbose.verbose(0)

# Go solve (transient)
md = pyissm.model.execute.solve(md, "tr")

# ---------------------------------------------------------------------
# Fields and tolerances to track changes (MATLAB-style regression checks)
# ---------------------------------------------------------------------
field_names = ["Gradient", "Misfit", "Rheology"]
field_tolerances = [2e-12, 1e-12, 1e-12]

# Typically pyISSM will store these in the transient solution results.
# Gradient index naming can differ (Gradient1 vs Gradient0, etc).
sol0 = md.results.TransientSolution[0]

field_values = [
    getattr(sol0, "Gradient1", None),               # gradient for first control
    getattr(sol0, "J", None),                       # total cost
    getattr(sol0, "MaterialsRheologyBbar", None),    # recovered Bbar time series
]

print("field_names:", field_names)
print("field_tolerances:", field_tolerances)
print("field_values types:", [type(v) for v in field_values])
