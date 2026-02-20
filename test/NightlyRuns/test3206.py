# Test Name: SquareShelfTransientCalibrationWithParamcodipack
import numpy as np
import pyissm

# -----------------------------
# 1) Generate observations (forward "truth" run)
# -----------------------------
md = pyissm.model.Model()
md = pyissm.model.mesh.triangle(md, "../assets/Exp/Square.exp", 50000.0)
md = pyissm.model.param.set_mask(md, "all", None)
md = pyissm.model.param.parameterize(md, "../assets/Par/SquareShelf.py")
md = pyissm.model.param.set_flow_equation(md, SSA="all")
md.cluster.np = 2

# Create real time series for B (vertex-based)
md.timestepping.interp_forcing = 0
md.timestepping.final_time = 2.0 * md.timestepping.time_step

nv = md.mesh.numberofvertices
B = 1.8e8 * np.ones((nv, 2))
B[np.where(md.mesh.x < md.mesh.y)[0], 1] = 1.4e8
md.materials.rheology_B = np.vstack([B, np.array([0.01, 2.0 * md.timestepping.time_step])])

# Initial values
md.initialization.vx = np.zeros(nv)
md.initialization.vy = np.zeros(nv)
md.initialization.pressure = np.zeros(nv)
md.initialization.temperature = np.zeros(nv)
md.basalforcings.geothermalflux = np.zeros(nv)
md.thermal.spctemperature = np.nan * np.ones(nv)

# Param basal forcings
md.basalforcings = pyissm.model.classes.basalforcings.linear()
md.basalforcings.deepwater_melting_rate = 50.0
md.basalforcings.deepwater_elevation = -500.0
md.basalforcings.upperwater_melting_rate = 0.0
md.basalforcings.upperwater_elevation = 0.0
md.basalforcings.groundedice_melting_rate = np.zeros(nv)
md.basalforcings.perturbation_melting_rate = np.zeros(nv)

md.transient.isthermal = 0

md = pyissm.model.execute.solve(md, "tr")

# -----------------------------
# 2) Set cost functions (one per transient time) + transient square misfits
# -----------------------------
from pyissm.model.classes.dependent import dependent
from pyissm.model.classes.independent import independent

# cost function classes (paths can differ; adjust if your repo uses different modules)
from pyissm.model.classes.cfsurface import cfsurfacelogvel
from pyissm.model.classes.cfsurface import cfsurfacesquaretransient

# Ensure containers exist
if getattr(md.outputdefinition, "definitions", None) is None:
    md.outputdefinition.definitions = []
if getattr(md.autodiff, "dependents", None) is None:
    md.autodiff.dependents = []
if getattr(md.autodiff, "independents", None) is None:
    md.autodiff.independents = []

weights = np.ones(nv)

# Iterate transient steps robustly (pyISSM can expose .steps or be directly iterable)

count = 1
for i in range(0, len(md.results.TransientSolution.steps)):
    sol = md.results.TransientSolution[i]
    vx_obs = getattr(sol, "Vx")
    vy_obs = getattr(sol, "Vy")
    time = sol.time

    # IMPORTANT: pyISSM constructors usually don't accept kwargs -> set attributes after init
    cf = cfsurfacelogvel()
    cf.name = f"LogVelMis{count}"
    cf.definitionstring = f"Outputdefinition{count}"
    cf.vxobs_string = "VxObs"
    cf.vxobs = vx_obs
    cf.vyobs_string = "VyObs"
    cf.vyobs = vy_obs
    cf.weights = weights
    cf.weights_string = "WeightsSurfaceObservation"
    cf.datatime = time

    md.outputdefinition.definitions.append(cf)

    dep = dependent()
    dep.name = f"Outputdefinition{count}"
    dep.type = "scalar"
    dep.fos_reverse_index = 1
    md.autodiff.dependents.append(dep)

    count += 1

# --- Deal with Vx separately (transient square misfit) ---
# MATLAB: vx_obs  = [[TransientSolution(:).Vx]/yts; [TransientSolution(:).time]];
vx_list = []
vy_list = []
surf_list = []
t_list = []
for i in range(0, len(md.results.TransientSolution.steps)):
    sol = md.results.TransientSolution[i]
    vx_list.append(getattr(sol, "Vx") / md.constants.yts)
    vy_list.append(getattr(sol, "Vy") / md.constants.yts)
    surf_list.append(getattr(sol, "Surface"))
    t = getattr(sol, "time", None)
    t_list.append(t)

vx_obs = np.vstack([np.column_stack(vx_list), np.asarray(t_list)])
vy_obs = np.vstack([np.column_stack(vy_list), np.asarray(t_list)])
surf_obs = np.vstack([np.column_stack(surf_list), np.asarray(t_list)])

# weights: [ones(nv,1); 0]
w_ts = np.concatenate([np.ones(nv), np.array([0.0])])

cfvx = cfsurfacesquaretransient()
cfvx.name = "VxMisfit_Transient"
cfvx.definitionstring = f"Outputdefinition{count}"
cfvx.model_string = "Vx"
cfvx.observations_string = "VxObs"
cfvx.observations = vx_obs
cfvx.weights = 500.0 * w_ts
cfvx.weights_string = "WeightsSurfaceObservation"
md.outputdefinition.definitions.append(cfvx)

dep = dependent()
dep.name = f"Outputdefinition{count}"
dep.type = "scalar"
dep.fos_reverse_index = 1
md.autodiff.dependents.append(dep)
count += 1

cfvy = cfsurfacesquaretransient()
cfvy.name = "VyMisfit_Transient"
cfvy.definitionstring = f"Outputdefinition{count}"
cfvy.model_string = "Vy"
cfvy.observations_string = "VyObs"
cfvy.observations = vy_obs
cfvy.weights = w_ts
cfvy.weights_string = "WeightsSurfaceObservation"
md.outputdefinition.definitions.append(cfvy)

dep = dependent()
dep.name = f"Outputdefinition{count}"
dep.type = "scalar"
dep.fos_reverse_index = 1
md.autodiff.dependents.append(dep)
count += 1

cfs = cfsurfacesquaretransient()
cfs.name = "SurfMisfit_Transient"
cfs.definitionstring = f"Outputdefinition{count}"
cfs.model_string = "Surface"
cfs.observations_string = "SurfaceObservation"
cfs.observations = surf_obs
cfs.weights = w_ts / md.constants.yts
cfs.weights_string = "WeightsSurfaceObservation"
md.outputdefinition.definitions.append(cfs)

dep = dependent()
dep.name = f"Outputdefinition{count}"
dep.type = "scalar"
dep.fos_reverse_index = 1
dep.nods = md.mesh.numberofvertices
md.autodiff.dependents.append(dep)
count += 1

# -----------------------------
# 3) Independents (Bbar + deepwater melt rate)
# -----------------------------
# Make B constant (except last time row)
md.materials.rheology_B[:-1, :] = 1.8e8

# Bounds for Bbar via cuffey(T)
from pyissm.tools.materials import cuffey

min_params = md.materials.rheology_B.copy()
max_params = md.materials.rheology_B.copy()
min_params[:-1, :] = cuffey(273)
max_params[:-1, :] = cuffey(200)

ind1 = independent()
ind1.name = "MaterialsRheologyBbar"
ind1.control_size = md.materials.rheology_B.shape[1]
ind1.type = "vertex"
ind1.min_parameters = min_params
ind1.max_parameters = max_params
ind1.control_scaling_factor = 1e8
md.autodiff.independents.append(ind1)

# Deepwater melting rate control
md.basalforcings.deepwater_melting_rate = 1.0
field = md.basalforcings.deepwater_melting_rate / md.constants.yts
name = "BasalforcingsDeepwaterMeltingRate"
scaling = 50.0 / md.constants.yts

ind2 = independent()
ind2.name = name
ind2.type = "vertex"
ind2.nods = md.mesh.numberofvertices
# MATLAB uses size(field,2); here scalar -> 1
ind2.control_size = 1
ind2.min_parameters = 1e-5 * field
ind2.max_parameters = 100.0 * field
ind2.control_scaling_factor = scaling
md.autodiff.independents.append(ind2)

# -----------------------------
# 4) Inversion / autodiff settings
# -----------------------------
from pyissm.model.classes.adm1qn3inversion import adm1qn3

md.inversion = adm1qn3(md.inversion)
md.inversion.iscontrol = 1
md.inversion.maxiter = 3
md.inversion.maxsteps = md.inversion.maxiter
md.inversion.dxmin = 1e-5

md.autodiff.isautodiff = 1
md.autodiff.driver = "fos_reverse"

md.settings.checkpoint_frequency = 2

# -----------------------------
# 5) Go solve (control run)
# -----------------------------
md = pyissm.model.execute.solve(md, "tr")

# -----------------------------
# 6) Fields and tolerances to track changes
# -----------------------------
# First transient step
sol0 = md.results.TransientSolution
if getattr(sol0, "steps", None) is not None:
    sol0 = sol0.steps[0]
else:
    sol0 = sol0[0]

field_names = ["Gradient1", "Gradient2", "Misfit", "Rheology", "DeepMelt"]
field_tolerances = [1e-10, 1e-10, 1e-10, 1e-10, 1e-10]
field_values = [
    getattr(sol0, "Gradient1"),
    getattr(sol0, "Gradient2"),
    getattr(sol0, "J"),
    getattr(sol0, "MaterialsRheologyBbar"),
    getattr(sol0, "BasalforcingsDeepwaterMeltingRate"),
]