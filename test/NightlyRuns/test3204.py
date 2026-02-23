# Test Name: SquareShelfTransientCalibrationNBEcodipackCheckpoint
# Same as test3201 but with checkpointing = 1

import numpy as np
import pyissm

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
md = pyissm.model.param.set_mask(md, "all", None)
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


count = 1

# Loop over transient results and add cost functions + dependents
for i in range(0, len(md.results.TransientSolution.steps)):
    sol = md.results.TransientSolution[i]
    vx_obs = sol.Vx
    vy_obs = sol.Vy
    z_obs = sol.Surface
    t_obs = sol.time

    weights = np.ones((nv, 1))

    # 1) Log velocity misfit
    cf = pyissm.model.classes.cf.cfsurfacelogvel()
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

    dep = pyissm.model.classes.dependent()
    dep.name = cf.definitionstring
    dep.type = "scalar"
    dep.fos_reverse_index = 1
    dep.nods = md.mesh.numberofvertices
    md.autodiff.dependents.append(dep)
    count += 1

    # 2) Vy square misfit (note MATLAB divides by yts for Vx/Vy observations)
    cf = pyissm.model.classes.cf.cfsurfacesquare()
    cf.name = f"VyMisfit{count+1}"
    cf.definitionstring = f"Outputdefinition{count+1}"
    cf.model_string = "Vy"
    cf.observation_string = "VyObs"
    cf.observation = vy_obs / md.constants.yts
    cf.weights = weights
    cf.weights_string = "WeightsSurfaceObservation"
    cf.datatime = t_obs
    md.outputdefinition.definitions.append(cf)

    dep = pyissm.model.classes.dependent()
    dep.name = cf.definitionstring
    dep.type = "scalar"
    dep.fos_reverse_index = 1
    dep.nods = md.mesh.numberofvertices
    md.autodiff.dependents.append(dep)
    count += 1

    # 3) Vx square misfit (MATLAB uses 500*weights here)
    cf = pyissm.model.classes.cf.cfsurfacesquare()
    cf.name = f"VxMisfit{count+1}"
    cf.definitionstring = f"Outputdefinition{count+1}"
    cf.model_string = "Vx"
    cf.observation_string = "VxObs"
    cf.observation = vx_obs / md.constants.yts
    cf.weights = 500.0 * weights
    cf.weights_string = "WeightsSurfaceObservation"
    cf.datatime = t_obs
    md.outputdefinition.definitions.append(cf)

    dep = pyissm.model.classes.dependent()
    dep.name = cf.definitionstring
    dep.type = "scalar"
    dep.fos_reverse_index = 1
    dep.nods = md.mesh.numberofvertices
    md.autodiff.dependents.append(dep)
    count += 1

    # 4) DEM / surface square misfit (MATLAB uses (1/yts)*weights)
    cf = pyissm.model.classes.cf.cfsurfacesquare()
    cf.name = f"DEMMisfit{count+1}"
    cf.definitionstring = f"Outputdefinition{count+1}"
    cf.model_string = "Surface"
    cf.observation_string = "SurfaceObservation"
    cf.observation = z_obs
    cf.weights = (1.0 / md.constants.yts) * weights
    cf.weights_string = "WeightsSurfaceObservation"
    cf.datatime = t_obs
    md.outputdefinition.definitions.append(cf)

    dep = pyissm.model.classes.dependent()
    dep.name = cf.definitionstring
    dep.type = "scalar"
    dep.fos_reverse_index = 1
    dep.nods = md.mesh.numberofvertices
    md.autodiff.dependents.append(dep)
    count += 1

# Independent (control)
min_params = md.materials.rheology_B.copy()
max_params = md.materials.rheology_B.copy()

# MATLAB: min_params(1:end-1,:) = cuffey(273); max_params(1:end-1,:) = cuffey(200);
min_params[:-1, :] = pyissm.tools.materials.cuffey(273.0)
max_params[:-1, :] = pyissm.tools.materials.cuffey(200.0)

ind = pyissm.model.classes.independent()
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
md.inversion = pyissm.model.classes.inversion.adm1qn3(other=md.inversion)
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
