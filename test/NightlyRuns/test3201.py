# Test Name: SquareShelfTransientCalibrationNBEcodipack
import numpy as np
import pyissm

# -----------------------------
# 1) Generate observations
# -----------------------------
md = pyissm.model.Model()
md = pyissm.model.mesh.triangle(md, "../assets/Exp/Square.exp", 100000.0)
md = pyissm.model.param.set_mask(md, "all", None)          # MATLAB: setmask(md,'all','')
md = pyissm.model.param.parameterize(md, "../assets/Par/SquareShelf.py")
md = pyissm.model.param.set_flow_equation(md, SSA="all")  # MATLAB: setflowequation(md,'SSA','all')

md.cluster.np = 2

# Create real time series for B
md.timestepping.interp_forcing = 0
md.timestepping.final_time = 2.0 * md.timestepping.time_step

ne = md.mesh.numberofelements
md.autodiff.isautodiff = False
# element-center coords for the "xbar < ybar" region
# MATLAB: mean(md.mesh.x(md.mesh.elements),2) < mean(md.mesh.y(md.mesh.elements),2)
elem = md.mesh.elements.astype(int) - 1  # ISSM elements are usually 1-based in MATLAB
xbar = md.mesh.x[elem].mean(axis=1)
ybar = md.mesh.y[elem].mean(axis=1)
region = xbar < ybar

B = 1.8e8 * np.ones((ne, 2))
B[region, 1] = 1.4e8

# MATLAB appends a "time row" for transient forcing:
# md.materials.rheology_B=[md.materials.rheology_B; 0.01 2*dt];
dt = md.timestepping.time_step
B = np.vstack([B, np.array([0.01, 2.0 * dt])])

md.materials.rheology_B = B

# Initial values
nv = md.mesh.numberofvertices
md.initialization.vx = np.zeros(nv)
md.initialization.vy = np.zeros(nv)
md.initialization.pressure = np.zeros(nv)
md.initialization.temperature = np.zeros(nv)

md.basalforcings.geothermalflux = np.zeros(nv)
md.thermal.spctemperature = np.full(nv, np.nan)

# md.toolkits.DefaultAnalysis = {"toolkit": "issm"}
md.toolkits.DefaultAnalysis = {
    "toolkit": "issm",
    "mat_type": "mpisparse",
    "vec_type": "mpi",
    "solver_type": "mumps"
}

# Solve once to generate the synthetic observations
print('toolkits:', md.toolkits.DefaultAnalysis, flush=True)
md = pyissm.model.execute.solve(md, "tr")

# -----------------------------
# 2) Modify rheology: now constant
# -----------------------------
md.materials.rheology_B[:-1, :] = 1.8e8

# -----------------------------
# 3) Set cost functions (one set per transient time)
# -----------------------------
from pyissm.model.classes.dependent import dependent
from pyissm.model.classes.independent import independent
from pyissm.model.classes.cfsurface import cfsurfacelogvel, cfsurfacesquare

def set_if(obj, **kvs):
    """Set attributes if they exist on this build's wrapper."""
    for k, v in kvs.items():
        if hasattr(obj, k):
            setattr(obj, k, v)

# Ensure containers exist
if getattr(md.outputdefinition, "definitions", None) is None:
    md.outputdefinition.definitions = []
if getattr(md.autodiff, "dependents", None) is None:
    md.autodiff.dependents = []

count = 0
for i in range(1, len(md.results.TransientSolution.steps)):
    vx_obs = md.results.TransientSolution[i].Vx
    vy_obs = md.results.TransientSolution[i].Vy
    z_obs  = md.results.TransientSolution[i].Surface
    t_obs  = md.results.TransientSolution[i].time

    weights = np.ones(nv)

    # ---- LogVel misfit ----
    cf = cfsurfacelogvel()
    set_if(cf,
        name=f"LogVelMis{count}",
        definitionstring=f"Outputdefinition{count}",
        vxobs_string="VxObs", vxobs=vx_obs,
        vyobs_string="VyObs", vyobs=vy_obs,
        weights=weights,
        weights_string="WeightsSurfaceObservation",
        datatime=t_obs,
    )
    md.outputdefinition.definitions.append(cf)

    dep = dependent()
    set_if(dep, name=f"Outputdefinition{count}", type="scalar", fos_reverse_index=1, nods=md.mesh.numberofvertices)
    md.autodiff.dependents.append(dep)
    count += 1

    # ---- Vy square misfit ----
    cf = cfsurfacesquare()
    set_if(cf,
        name=f"VyMisfit{count}",
        definitionstring=f"Outputdefinition{count}",
        model_string="Vy",
        observation_string="VyObs",
        observation=vy_obs / md.constants.yts,
        weights=weights,
        weights_string="WeightsSurfaceObservation",
        datatime=t_obs,
    )
    md.outputdefinition.definitions.append(cf)

    dep = dependent()
    set_if(dep, name=f"Outputdefinition{count}", type="scalar", fos_reverse_index=1,         nods = md.mesh.numberofvertices,)
    md.autodiff.dependents.append(dep)
    count += 1

    # ---- Vx square misfit (500*weights) ----
    cf = cfsurfacesquare()
    set_if(cf,
        name=f"VxMisfit{count}",
        definitionstring=f"Outputdefinition{count}",
        model_string="Vx",
        observation_string="VxObs",
        observation=vx_obs / md.constants.yts,
        weights=500.0 * weights,
        weights_string="WeightsSurfaceObservation",
        datatime=t_obs,
    )
    md.outputdefinition.definitions.append(cf)

    dep = dependent()
    set_if(dep, name=f"Outputdefinition{count}", type="scalar", fos_reverse_index=1, nods=md.mesh.numberofvertices)
    md.autodiff.dependents.append(dep)
    count += 1

    # ---- DEM / Surface misfit (1/yts * weights) ----
    cf = cfsurfacesquare()
    set_if(cf,
        name=f"DEMMisfit{count}",
        definitionstring=f"Outputdefinition{count}",
        model_string="Surface",
        observation_string="SurfaceObservation",
        observation=z_obs,
        weights=(1.0 / md.constants.yts) * weights,
        weights_string="WeightsSurfaceObservation",
        datatime=t_obs,  # test that this optional argument is accepted by the wrapper
    )
    md.outputdefinition.definitions.append(cf)

    dep = dependent()
    set_if(dep, name=f"Outputdefinition{count}", type="scalar", fos_reverse_index=1,nods = md.mesh.numberofvertices)
    md.autodiff.dependents.append(dep)
    count += 1


# -----------------------------
# 4) Independent: MaterialsRheologyBbar
# -----------------------------
from pyissm.tools.materials import cuffey

min_params = md.materials.rheology_B.copy()
max_params = md.materials.rheology_B.copy()
min_params[:-1, :] = cuffey(273)
max_params[:-1, :] = cuffey(200)

md.autodiff.independents = []

ind = independent()
set_if(ind,
    name="MaterialsRheologyBbar",
    control_size=md.materials.rheology_B.shape[1],
    type="vertex",  # if your build wants element-based control, switch to "element"
    min_parameters=min_params,
    max_parameters=max_params,
    control_scaling_factor=1e8,
)
md.autodiff.independents = [ind]

# -----------------------------
# 5) Inversion / AD settings and solve
# -----------------------------
# MATLAB: md.inversion=adm1qn3inversion(md.inversion);

md.inversion = pyissm.model.classes.inversion.adm1qn3(md.inversion)


md.inversion.iscontrol = 1
md.inversion.maxiter = 3
md.inversion.maxsteps = md.inversion.maxiter
md.inversion.dxmin = 1e-5

md.autodiff.isautodiff = True
md.autodiff.driver = "fos_reverse"


md = pyissm.model.execute.solve(md, "tr")

# -----------------------------
# 6) Fields and tolerances to track changes (like the MATLAB test harness)
# -----------------------------
field_names = ["Gradient", "Misfit", "Rheology"]
field_tolerances = [1e-12, 1e-12, 1e-12]
field_values = [
    md.results.TransientSolution[0].Gradient1,
    md.results.TransientSolution[0].J,
    md.results.TransientSolution[0].MaterialsRheologyBbar,
]
