# Test Name: SquareShelfTransientLevelsetMisfitcodipack
import os
import numpy as np
import pyissm

# -----------------------------
# 1) Mesh / mask / params / flow eqn / cluster
# -----------------------------
md = pyissm.model.Model()
md = pyissm.model.mesh.triangle(md, "../assets/Exp/Square.exp", 50000.0)
md = pyissm.model.param.set_mask(md, "all", None)
md = pyissm.model.param.parameterize(md, "../assets/Par/SquareShelf.py")
md = pyissm.model.param.set_flow_equation(md, SSA="all")

# MATLAB: md.cluster=generic('name',oshostname(),'np',3);
# pyISSM: usually just set number of MPI processes
md.cluster.np = 3

# -----------------------------
# 2) Levelset + timestep settings
# -----------------------------
# Do not kill icebergs as all is floating
md.levelset.kill_icebergs = 0

x = md.mesh.x
xmin = np.min(x)
xmax = np.max(x)
Lx = (xmax - xmin)
alpha = 2.0 / 3.0

# MATLAB: ((x - alpha*Lx)>0) - ((x - alpha*Lx)<0)
md.mask.ice_levelset = ((x - alpha * Lx) > 0).astype(float) - ((x - alpha * Lx) < 0).astype(float)

md.timestepping.time_step = 10.0
md.timestepping.final_time = 30.0

# Transient flags
md.transient.isstressbalance = 1
md.transient.ismasstransport = 1
md.transient.issmb = 1
md.transient.isthermal = 0
md.transient.isgroundingline = 0
md.transient.ismovingfront = 1

# -----------------------------
# 3) Calving / frontal forcings / levelset constraints
# -----------------------------
md.calving = pyissm.model.classes.calving.levermann()

md.calving.coeff = 4.89e13 * np.ones(md.mesh.numberofvertices)
md.frontalforcings.meltingrate = np.zeros(md.mesh.numberofvertices)

md.levelset.spclevelset = np.nan * np.ones(md.mesh.numberofvertices)
md.levelset.migration_max = 1e8

# -----------------------------
# 4) First transient solve ("truth run")
# -----------------------------
md = pyissm.model.execute.solve(md, "tr")

# -----------------------------
# 5) Modify rheology, now constant
# -----------------------------
# NOTE: in MATLAB you do (1:end-1,:) because the last row can be [time, ...] for time series forcing.
# Keep that same convention here.
md.materials.rheology_B[:-1] = 1.8e8

# -----------------------------
# 6) Set cost functions: one cflevelsetmisfit + one dependent per transient time
# -----------------------------
weights = np.ones(md.mesh.numberofvertices)

# Ensure containers exist (pyISSM sometimes uses lists, sometimes special containers)
if getattr(md.outputdefinition, "definitions", None) is None:
    md.outputdefinition.definitions = []
if getattr(md.autodiff, "dependents", None) is None:
    md.autodiff.dependents = []
if getattr(md.autodiff, "independents", None) is None:
    md.autodiff.independents = []

# Helpers: iterate transient steps robustly (pyISSM can expose .steps or be directly iterable)
ts = md.results.TransientSolution
steps = getattr(ts, "steps", None)
if steps is None:
    # some builds allow: for sol in md.results.TransientSolution:
    steps = ts

# Import classes (paths can differ across pyISSM builds)
from pyissm.model.classes.dependent import dependent
from pyissm.model.classes.independent import independent
from pyissm.model.classes.cflevelsetmisfit import cflevelsetmisfit

# reinitializelevelset function location can vary; try common places
from pyissm.model.classes.levelset import reinitialize

count = 1
for sol in steps:
    # Time: in some builds it’s sol.time, in others sol.Time
    time = getattr(sol, "time", None)
    if time is None:
        time = getattr(sol, "Time", None)

    # Modelled levelset at that step
    phi = getattr(sol, "MaskIceLevelset", None)
    if phi is None:
        raise AttributeError("TransientSolution step has no 'MaskIceLevelset' field.")

    obs = reinitialize(md, phi)

    # IMPORTANT: many pyISSM class constructors DO NOT accept kwargs like MATLAB.
    # Pattern that usually works: instantiate with no args, then set attributes.
    mis = cflevelsetmisfit()
    mis.name = f"LevelsetMisfit{count}"
    mis.definitionstring = f"Outputdefinition{count}"
    mis.model_string = "MaskIceLevelset"
    mis.observation_string = "LevelsetObservation"
    mis.observation = obs
    mis.weights = weights
    mis.weights_string = "WeightsLevelsetObservation"
    mis.datatime = time

    md.outputdefinition.definitions.append(mis)

    dep = dependent()
    dep.name = f"Outputdefinition{count}"
    dep.type = "scalar"
    dep.fos_reverse_index = 1
    dep.nods = md.mesh.numberofvertices
    md.autodiff.dependents.append(dep)

    count += 1

# -----------------------------
# 7) Independent (control) on MaterialsRheologyBbar with bounds from cuffey(T)
# -----------------------------
min_params = md.materials.rheology_B.copy()
max_params = md.materials.rheology_B.copy()

# cuffey(T) helper: location varies; try common import

from pyissm.model.materials import cuffey


if cuffey is None:
    raise ImportError("Could not import cuffey(T). Locate it in your pyISSM install and update the import.")

min_params[:-1, :] = cuffey(273)  # warm -> lower viscosity (check your cuffey convention)
max_params[:-1, :] = cuffey(200)  # cold -> higher viscosity

ind = independent()
ind.name = "MaterialsRheologyBbar"
ind.control_size = md.materials.rheology_B.shape[1]
ind.type = "vertex"  # matches your MATLAB comment
ind.min_parameters = min_params
ind.max_parameters = max_params
ind.control_scaling_factor = 1e8
md.autodiff.independents.append(ind)

# -----------------------------
# 8) Inversion + autodiff settings
# -----------------------------
# Ensure inversion object exists; then wrap/convert to adm1qn3inversion like MATLAB

from pyissm.model.classes.inversion import adm1qn3inversion  # alternate path some builds use

md.inversion = adm1qn3inversion(md.inversion)
md.inversion.iscontrol = 1
md.inversion.maxiter = 4
md.inversion.maxsteps = md.inversion.maxiter
md.inversion.dxmin = 1e-5

md.autodiff.isautodiff = 1
md.autodiff.driver = "fos_reverse"

# -----------------------------
# 9) Solve (control run)
# -----------------------------
md = pyissm.model.execute.solve(md, "tr")

# -----------------------------
# 10) Fields / tolerances (as in MATLAB end block)
# -----------------------------
field_names = ["Gradient", "Misfit", "Rheology"]
field_tolerances = [1e-12, 1e-12, 1e-12]

# first transient step result access (pyISSM sometimes uses [0] not (1))
sol0 = getattr(md.results, "TransientSolution", None)
if getattr(sol0, "steps", None) is not None:
    sol0 = sol0.steps[0]
else:
    sol0 = sol0[0]

field_values = [
    getattr(sol0, "Gradient1"),
    getattr(sol0, "J"),
    getattr(sol0, "MaterialsRheologyBbar"),
]