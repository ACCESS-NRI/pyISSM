# Test Name: SquareShelfTransientCalibrationNBEcodipackCheckpoint
# Same as test3201 but with checkpointing = 1

import numpy as np
import pyissm

# These imports may vary slightly depending on your tree layout
from pyissm.model import solve, mesh
from pyissm.model.param import set_mask, parameterize, set_flow_equation
from pyissm.model.classes.autodiff import dependent, independent
from pyissm.model.classes.outputdefinition import outputdefinition
from pyissm.model.classes.verbose import verbose

# Cost-function definition classes (assumed in your tree)
from pyissm.model.classes.cfsurfacelogvel import cfsurfacelogvel
from pyissm.model.classes.cfsurfacesquare import cfsurfacesquare

# Inversion driver (your cleaned pyISSM-style class)
from pyissm.model.classes.inversion import adm1qn3

# Temperature-to-B bounds helper (name/location may differ in your tree)
# If this import fails in your repo, replace with the correct path to cuffey().
from pyissm.model.materials import cuffey


def _element_centroids(md):
    """
    Return element centroid x,y arrays.

    Handles the common ISSM convention where md.mesh.elements can be 1-based.
    """
    elems = md.mesh.elements
    # If elements are 1-based (typical in ISSM), convert to 0-based for numpy indexing.
    if elems.min() == 1:
        elems0 = elems - 1
    else:
        elems0 = elems

    x = md.mesh.x[elems0].mean(axis=1)
    y = md.mesh.y[elems0].mean(axis=1)
    return x, y


# -------------------------
# Generate observations run
# -------------------------
md = pyissm.model.Model()
md = pyissm.model.mesh.triangle(md, "../assets/Exp/Square.exp", 100000.0)
md = pyissm.param.parameterize.set_mask(md, "all", None)
md = pyissm.model.param.parameterize(md, "../assets/Par/SquareShelf.py")
md = pyissm.model.param.set_flow_equation(md, SSA="all")
md.cluster.np = 2

# Create real time series for B
md.timestepping.interp_forcing = 0
md.timestepping.final_time = 2 * md.timestepping.time_step

ne = md.mesh.numberofelements
dt = md.timestepping.time_step

# MATLAB: md.materials.rheology_B = 1.8e8 * ones(ne,2);
md.materials.rheology_B = 1.8e8 * np.ones((ne, 2))

# MATLAB condition uses element centroid means:
# mean(x(elements),2) < mean(y(elements),2) -> set column 2 (MATLAB indexing) to 1.4e8
cx, cy = _element_centroids(md)
mask = cx < cy
md.materials.rheology_B[mask, 1] = 1.4e8  # python col=1 is MATLAB col=2

# MATLAB appends a time row: [0.01 2*dt]
md.materials.rheology_B = np.vstack([md.materials.rheology_B, np.array([0.01, 2 * dt])])

# Initial values
nv = md.mesh.numberofvertices
md.initialization.vx = np.zeros((nv, 1))
md.initialization.vy = np.zeros((nv, 1))
md.initialization.pressure = np.zeros((nv, 1))
md.initialization.temperature = np.zeros((nv, 1))
md.basalforcings.geothermalflux = np.zeros((nv, 1))
md.thermal.spctemperature = np.full((nv, 1), np.nan)

# Solve transient to generate "observations"
md = pyissm.model.execute.solve(md, "tr")

# -------------------------
# Inversion setup run
# -------------------------

# Modify rheology, now constant (keep the final appended time row)
md.materials.rheology_B[:-1, :] = 1.8e8

# Make sure outputdefinition/autodiff containers exist
if not hasattr(md, "outputdefinition") or md.outputdefinition is None:
    md.outputdefinition = outputdefinition()
if not hasattr(md.outputdefinition, "definitions") or md.outputdefinition.definitions is None:
    md.outputdefinition.definitions = []

if not hasattr(md, "autodiff") or md.autodiff is None:
    md.autodiff = pyissm.model.classes.autodiff.autodiff()
if not hasattr(md.autodiff, "dependents") or md.autodiff.dependents is None:
    md.autodiff.dependents = []
if not hasattr(md.autodiff, "independents") or md.autodiff.independents is None:
    md.autodiff.independents = []

count = 0

# Loop over transient results and add cost functions + dependents
for sol in md.results.TransientSolution:
    vx_obs = sol.Vx
    vy_obs = sol.Vy
    z_obs = sol.Surface
    t_obs = sol.time

    weights = np.ones((nv, 1))

    # 1) Log velocity misfit
    cf = cfsurfacelogvel()
    cf.name = f"LogVelMis{count+1}"
    cf.definitionstring = f"Outputdefinition{count+1}"
    cf.vxobs_string = "VxObs"
    cf.vxobs = vx_obs
    cf.vyobs_string = "VyObs"
    cf.vyobs = vy_obs
    cf.weights = weights
    cf.weights_string = "WeightsSurfaceObservation"
    cf.datatime = t_obs
    md.outputdefinition.definitions.append(cf)

    dep = dependent()
    dep.name = cf.definitionstring
    dep.type = "scalar"
    dep.fos_reverse_index = 1
    md.autodiff.dependents.append(dep)
    count += 1

    # 2) Vy square misfit (note MATLAB divides by yts for Vx/Vy observations)
    cf = cfsurfacesquare()
    cf.name = f"VyMisfit{count+1}"
    cf.definitionstring = f"Outputdefinition{count+1}"
    cf.model_string = "Vy"
    cf.observation_string = "VyObs"
    cf.observation = vy_obs / md.constants.yts
    cf.weights = weights
    cf.weights_string = "WeightsSurfaceObservation"
    cf.datatime = t_obs
    md.outputdefinition.definitions.append(cf)

    dep = dependent()
    dep.name = cf.definitionstring
    dep.type = "scalar"
    dep.fos_reverse_index = 1
    md.autodiff.dependents.append(dep)
    count += 1

    # 3) Vx square misfit (MATLAB uses 500*weights here)
    cf = cfsurfacesquare()
    cf.name = f"VxMisfit{count+1}"
    cf.definitionstring = f"Outputdefinition{count+1}"
    cf.model_string = "Vx"
    cf.observation_string = "VxObs"
    cf.observation = vx_obs / md.constants.yts
    cf.weights = 500.0 * weights
    cf.weights_string = "WeightsSurfaceObservation"
    cf.datatime = t_obs
    md.outputdefinition.definitions.append(cf)

    dep = dependent()
    dep.name = cf.definitionstring
    dep.type = "scalar"
    dep.fos_reverse_index = 1
    md.autodiff.dependents.append(dep)
    count += 1

    # 4) DEM / surface square misfit (MATLAB uses (1/yts)*weights)
    cf = cfsurfacesquare()
    cf.name = f"DEMMisfit{count+1}"
    cf.definitionstring = f"Outputdefinition{count+1}"
    cf.model_string = "Surface"
    cf.observation_string = "SurfaceObservation"
    cf.observation = z_obs
    cf.weights = (1.0 / md.constants.yts) * weights
    cf.weights_string = "WeightsSurfaceObservation"
    cf.datatime = t_obs
    md.outputdefinition.definitions.append(cf)

    dep = dependent()
    dep.name = cf.definitionstring
    dep.type = "scalar"
    dep.fos_reverse_index = 1
    md.autodiff.dependents.append(dep)
    count += 1

# Independent (control)
min_params = md.materials.rheology_B.copy()
max_params = md.materials.rheology_B.copy()

# MATLAB: min_params(1:end-1,:) = cuffey(273); max_params(1:end-1,:) = cuffey(200);
min_params[:-1, :] = cuffey(273.0)
max_params[:-1, :] = cuffey(200.0)

ind = independent()
ind.name = "MaterialsRheologyBbar"
ind.control_size = md.materials.rheology_B.shape[1]
ind.type = "vertex"  # MATLAB comment says “Really needed??” — keep parity
ind.min_parameters = min_params
ind.max_parameters = max_params
ind.control_scaling_factor = 1e8

# Put it in the first slot
if len(md.autodiff.independents) == 0:
    md.autodiff.independents.append(ind)
else:
    md.autodiff.independents[0] = ind

# Inversion driver: ADM1QN3
md.inversion = adm1qn3(other=md.inversion)
md.inversion.iscontrol = 1
md.inversion.maxiter = 3
md.inversion.maxsteps = md.inversion.maxiter
md.inversion.dxmin = 1e-5

# Autodiff / checkpointing
md.autodiff.isautodiff = True
md.autodiff.driver = "fos_reverse"
md.settings.checkpoint_frequency = 1

md = pyissm.model.execute.solve(md, "tr")

# Fields and tolerances to track changes
field_names = ["Gradient", "Misfit", "Rheology"]
field_tolerances = [1e-12, 1e-12, 1e-12]
field_values = [
    md.results.TransientSolution[0].Gradient1,
    md.results.TransientSolution[0].J,
    md.results.TransientSolution[0].MaterialsRheologyBbar,
]
