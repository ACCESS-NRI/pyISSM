# Test Name: SquareShelfTransientLevelsetMisfitcodipack
import numpy as np
import pyissm

# --- Mesh / param / flow eqn / cluster ---
md = pyissm.model.Model()
md = pyissm.model.mesh.triangle(pyissm.model.Model(), "../assets/Exp/Square.exp", 50000.0)
md = pyissm.model.param.set_mask(md, "all", None)
md = pyissm.model.param.parameterize(md, "../assets/Par/SquareShelf.py")
md = pyissm.model.param.set_flow_equation(md, SSA="all")

# MATLAB: md.cluster=generic('name',oshostname(),'np',3);
# In pyISSM you typically just set np (or use a cluster helper if you have one configured)
md.cluster.np = 3

# --- Levelset options ---
# Do not kill icebergs as all is floating
md.levelset.kill_icebergs = 0

x = md.mesh.x
xmin = np.min(x)
xmax = np.max(x)
Lx = (xmax - xmin)
alpha = 2.0 / 3.0

# MATLAB: ((x - alpha*Lx)>0) - ((x - alpha*Lx)<0);
# -> +1 on one side, -1 on the other, 0 exactly on the interface
phi = x - alpha * Lx
md.mask.ice_levelset = (phi > 0).astype(float) - (phi < 0).astype(float)
md.materials.rheology_B[md.mesh.x < md.mesh.y, 1] = 1.4e8
md.materials.rheology_B = np.vstack([md.materials.rheology_B, [0.01, 2*md.timestepping.time_step]])

# --- Time stepping ---
md.timestepping.time_step = 10.0
md.timestepping.final_time = 30.0

# --- Transient flags ---
md.transient.isstressbalance = 1
md.transient.ismasstransport = 1
md.transient.issmb = 1
md.transient.isthermal = 0
md.transient.isgroundingline = 0
md.transient.ismovingfront = 1

# --- Calving / frontal forcings / levelset constraints ---
md.calving = pyissm.model.classes.calving.levermann()
md.calving.coeff = 4.89e13 * np.ones(md.mesh.numberofvertices)

md.frontalforcings.meltingrate = np.zeros(md.mesh.numberofvertices)

md.levelset.spclevelset = np.full(md.mesh.numberofvertices, np.nan)
md.levelset.migration_max = 1e8

# --- Forward transient solve (truth run) ---
md = pyissm.model.execute.solve(md, "tr")

# --- Modify rheology, now constant ---
# Now modify all rows except the last (time) row
md.materials.rheology_B[:-1, :] = 1.8e8



  # keep the final time row


# --- Cost function setup over all transient times ---
weights = np.ones(md.mesh.numberofvertices,)

# Make sure containers exist
if md.outputdefinition.definitions is None:
    md.outputdefinition.definitions = []
if md.autodiff.dependents is None:
    md.autodiff.dependents = []

# reinitializelevelset(md, levelset) equivalent:
# (name can differ depending on your pyISSM version; keep one of these)
def _reinit_levelset(md, ls):
    # try a couple of common spellings/locations
    if hasattr(pyissm, "reinitializelevelset"):
        return pyissm.reinitializelevelset(md, ls)
    if hasattr(pyissm.model, "levelset") and hasattr(pyissm.model.levelset, "reinitializelevelset"):
        return pyissm.model.levelset.reinitializelevelset(md, ls)
    if hasattr(pyissm.model, "levelset") and hasattr(pyissm.model.levelset, "reinitialize_levelset"):
        return pyissm.model.levelset.reinitialize_levelset(md, ls)
    # fallback: use as-is (still lets you build the test)
    return ls

count = 1
for i in range(0, len(md.results.TransientSolution.steps)):
    sol = md.results.TransientSolution[i]
    time = sol.time

    obs = _reinit_levelset(md, sol.MaskIceLevelset)

    cf = pyissm.model.classes.cflevelsetmisfit()

    cf.name = f"LevelsetMisfit{count}"
    cf.definitionstring = f"Outputdefinition{count}"
    cf.model_string = "MaskIceLevelset"
    cf.observation_string = "LevelsetObservation"
    cf.observation = obs
    cf.weights = weights
    cf.weights_string = "WeightsLevelsetObservation"
    cf.datatime = time

    md.outputdefinition.definitions.append(cf)
    
    dep = pyissm.model.classes.dependent()
    dep.name = f"Outputdefinition{count}"
    dep.type = "scalar"
    dep.fos_reverse_index = 1
    dep.nods = md.mesh.numberofvertices
    md.autodiff.dependents.append(dep)

    count += 1

# --- Independent (control) bounds ---
# MATLAB:
#   min_params = md.materials.rheology_B; min_params(1:end-1,:) = cuffey(273);
#   max_params = md.materials.rheology_B; max_params(1:end-1,:) = cuffey(200);
min_params = md.materials.rheology_B.copy()
max_params = md.materials.rheology_B.copy()

# cuffey(T) should exist in your ISSM/pyISSM environment; otherwise replace with your own mapping.
# (In ISSM, cuffey() returns B(T) typically.)
min_params[:-1, :] = pyissm.tools.materials.cuffey(273)
max_params[:-1, :] = pyissm.tools.materials.cuffey(200)

if md.autodiff.independents is None:
    md.autodiff.independents = []

# --- Add independent (control) ---
ind = pyissm.model.classes.independent()
ind.name = "MaterialsRheologyBbar"
ind.control_size = md.materials.rheology_B.shape[1]
ind.type = "vertex"  
ind.min_parameters = min_params
ind.max_parameters = max_params
ind.control_scaling_factor = 1e8

md.autodiff.independents = [ind]

# --- Inversion / autodiff settings ---
md.inversion = pyissm.model.classes.inversion.adm1qn3(md.inversion)
md.inversion.iscontrol = 1
md.inversion.maxiter = 4
md.inversion.maxsteps = md.inversion.maxiter
md.inversion.dxmin = 1e-5

md.autodiff.isautodiff = True
md.autodiff.driver = "fos_reverse"

# --- Go solve (control run) ---
md = pyissm.model.execute.solve(md, "tr")

g = md.results.TransientSolution[0].Gradient1
if g.shape == (2, md.mesh.numberofvertices):
    g = g[0, :].reshape((-1, 1))      # take control 1

Bbar = md.results.TransientSolution[0].MaterialsRheologyBbar
if Bbar.shape == (2, md.mesh.numberofvertices):
    Bbar = Bbar[0, :].reshape((-1, 1))


# --- Fields and tolerances to track changes ---
field_names = ["Gradient", "Misfit", "Rheology"]
field_tolerances = [1e-12, 1e-12, 1e-12]
field_values = [
    g,
    md.results.TransientSolution[0].J,
    Bbar,
]
