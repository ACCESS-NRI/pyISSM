# Test Name: SquareShelfTransientCalibrationWithParamcodipack

import numpy as np
import pyissm

# -----------------------------
# Generate observations
# -----------------------------
md = pyissm.model.Model()
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 50000.)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelf.py')
md = pyissm.model.param.set_flow_equation(md, SSA='all')
md.cluster.np = 2

# -----------------------------
# Create real time series for B (vertex-based)
# -----------------------------
md.timestepping.interp_forcing = 0
md.timestepping.final_time = 2.0 * md.timestepping.time_step

nv = md.mesh.numberofvertices
md.materials.rheology_B = 1.8e8 * np.ones((md.mesh.numberofvertices, 2))
md.materials.rheology_B[md.mesh.x < md.mesh.y, 1] = 1.4e8
md.materials.rheology_B = np.vstack([md.materials.rheology_B, [0.01, 2*md.timestepping.time_step]])
# -----------------------------
# Initial values
# -----------------------------
md.initialization.vx = np.zeros((nv,))
md.initialization.vy = np.zeros((nv,))
md.initialization.pressure = np.zeros((nv,))
md.initialization.temperature = np.zeros((nv,))
md.basalforcings.geothermalflux = np.zeros((nv,))
md.thermal.spctemperature = np.full((nv,), np.nan)

# -----------------------------
# Param: linear basal forcings
# -----------------------------
md.basalforcings = pyissm.model.classes.basalforcings.linear()
md.basalforcings.deepwater_melting_rate = 50.0              # m/yr ice equivalent
md.basalforcings.deepwater_elevation = -500.0
md.basalforcings.upperwater_melting_rate = 0.0              # no melting for zb>=0
md.basalforcings.upperwater_elevation = 0.0                 # sea level
md.basalforcings.groundedice_melting_rate = np.zeros((nv,)) # no melting on grounded ice
md.basalforcings.perturbation_melting_rate = np.zeros((nv,)) # no perturbation melting
md.transient.isthermal = 0

# Forward solve to create observations
md = pyissm.model.execute.solve(md, 'tr')

# -----------------------------
# Set cost function: per-time LogVel misfits
# -----------------------------
md.outputdefinition.definitions = []
md.autodiff.dependents = []

count = 1
for i in range(0, len(md.results.TransientSolution.steps)):
    sol = md.results.TransientSolution[i]
    vx_obs = sol.Vx
    vy_obs = sol.Vy
    time = sol.time
    weights = np.ones((nv,))

    od = pyissm.model.classes.cfsurface.cfsurfacelogvel()
    od.name = f'LogVelMis{count}'
    od.definitionstring = f'Outputdefinition{count}'
    od.vxobs_string = 'VxObs'
    od.vxobs = vx_obs
    od.vyobs_string = 'VyObs'
    od.vyobs = vy_obs
    od.weights = weights
    od.weights_string = 'WeightsSurfaceObservation'
    od.datatime = time
    md.outputdefinition.definitions.append(od)

    dep = pyissm.model.classes.dependent()
    dep.name = f'Outputdefinition{count}'
    dep.type = 'scalar'
    dep.fos_reverse_index = 1
    dep.nods = md.mesh.numberofvertices
    md.autodiff.dependents.append(dep)

    count += 1

# -----------------------------
# Deal with Vx/Vy/Surface as transient observations (stack fields + times)
# MATLAB builds: [ [fields]/yts ; [times] ]
# We'll build arrays of shape (nv+1, nt) with last row = times.
# -----------------------------
times = times = np.array([sol.time for sol in md.results.TransientSolution.steps])
nt = len(times)

vx_stack = np.column_stack([sol.Vx for sol in md.results.TransientSolution.steps]) / md.constants.yts
vy_stack = np.column_stack([sol.Vy for sol in md.results.TransientSolution.steps]) / md.constants.yts
surf_stack = np.column_stack([sol.Surface for sol in md.results.TransientSolution.steps])

vx_obs_tr = np.vstack([vx_stack, times])
vy_obs_tr = np.vstack([vy_stack, times])
surf_obs_tr = np.vstack([surf_stack, times])

weights_tr = np.vstack([np.ones((nv, 1)), np.array([[0.0]])])  # last row weight for time row = 0

# Vx transient square misfit
od_vx = pyissm.model.classes.cfsurface.cfsurfacesquaretransient()
od_vx.name = 'VxMisfit_Transient'
od_vx.definitionstring = f'Outputdefinition{count}'
od_vx.model_string = 'Vx'
od_vx.observations_string = 'VxObs'
od_vx.observations = vx_obs_tr
od_vx.weights = 500.0 * weights_tr
od_vx.weights_string = 'WeightsSurfaceObservation'
md.outputdefinition.definitions.append(od_vx)

dep_vx = pyissm.model.classes.dependent()
dep_vx.name = f'Outputdefinition{count}'
dep_vx.type = 'scalar'
dep_vx.fos_reverse_index = 1
dep_vx.nods = md.mesh.numberofvertices
md.autodiff.dependents.append(dep_vx)
count += 1

# Vy transient square misfit
od_vy = pyissm.model.classes.cfsurface.cfsurfacesquaretransient()
od_vy.name = 'VyMisfit_Transient'
od_vy.definitionstring = f'Outputdefinition{count}'
od_vy.model_string = 'Vy'
od_vy.observations_string = 'VyObs'
od_vy.observations = vy_obs_tr
od_vy.weights = weights_tr
od_vy.weights_string = 'WeightsSurfaceObservation'
md.outputdefinition.definitions.append(od_vy)

dep_vy = pyissm.model.classes.dependent()
dep_vy.name = f'Outputdefinition{count}'
dep_vy.type = 'scalar'
dep_vy.fos_reverse_index = 1
dep_vy.nods = md.mesh.numberofvertices
md.autodiff.dependents.append(dep_vy)
count += 1

# Surface transient square misfit
od_surf = pyissm.model.classes.cfsurface.cfsurfacesquaretransient()
od_surf.name = 'SurfMisfit_Transient'
od_surf.definitionstring = f'Outputdefinition{count}'
od_surf.model_string = 'Surface'
od_surf.observations_string = 'SurfaceObservation'
od_surf.observations = surf_obs_tr
od_surf.weights = weights_tr / md.constants.yts
od_surf.weights_string = 'WeightsSurfaceObservation'
md.outputdefinition.definitions.append(od_surf)

dep_surf = pyissm.model.classes.dependent()
dep_surf.name = f'Outputdefinition{count}'
dep_surf.type = 'scalar'
dep_surf.fos_reverse_index = 1
dep_surf.nods = md.mesh.numberofvertices
md.autodiff.dependents.append(dep_surf)
count += 1

# -----------------------------
# Independents: (1) MaterialsRheologyBbar, (2) Deepwater melting rate
# -----------------------------
# reset B constant before inversion
# md.materials.rheology_B[:-1, :] = 1.8e8

min_params = md.materials.rheology_B.copy()
max_params = md.materials.rheology_B.copy()
min_params[:-1, :] = pyissm.tools.materials.cuffey(273)
max_params[:-1, :] = pyissm.tools.materials.cuffey(200)

indep1 = pyissm.model.classes.independent()
indep1.name = 'MaterialsRheologyBbar'
indep1.control_size = md.materials.rheology_B.shape[1]
indep1.type = 'vertex'
indep1.min_parameters = min_params
indep1.max_parameters = max_params
indep1.control_scaling_factor = 1e8

# second control: BasalforcingsDeepwaterMeltingRate (set deepwater rate to 1 m/yr first)
md.basalforcings.deepwater_melting_rate = 1.0  # m/yr ice equivalent
field = md.basalforcings.deepwater_melting_rate / md.constants.yts  # scalar (1/s)
name = 'BasalforcingsDeepwaterMeltingRate'
scaling = 50.0 / md.constants.yts

indep2 = pyissm.model.classes.independent()
indep2.name = name
indep2.type = 'vertex'
indep2.nods = md.mesh.numberofvertices
# control_size in MATLAB is size(field,2); for scalar -> 1
indep2.control_size = 1
indep2.min_parameters = 1e-5 * field
indep2.max_parameters = 100.0 * field
indep2.control_scaling_factor = scaling

md.autodiff.independents = [indep1, indep2]

# -----------------------------
# Inversion + autodiff + checkpointing
# -----------------------------
md.inversion = pyissm.model.classes.inversion.adm1qn3(md.inversion)
md.inversion.iscontrol = 1
md.inversion.maxiter = 3
md.inversion.maxsteps = md.inversion.maxiter
md.inversion.dxmin = 1e-5

md.autodiff.isautodiff = True
md.autodiff.driver = 'fos_reverse'

md.settings.checkpoint_frequency = 2

# Go solve!
md = pyissm.model.execute.solve(md, 'tr')

# -----------------------------
# Fields and tolerances to track changes
# -----------------------------
field_names = ['Gradient1', 'Gradient2', 'Misfit', 'Rheology', 'DeepMelt']
field_tolerances = [1e-10] * len(field_names)

ts1 = md.results.TransientSolution[0]  # MATLAB (1)
field_values = [
    ts1.Gradient1,
    ts1.Gradient2,
    ts1.J,
    ts1.MaterialsRheologyBbar,
    ts1.BasalforcingsDeepwaterMeltingRate
]

# -----------------------------
# Gradient validation block (MATLAB code after `return;`)
# Kept here as an optional function you can run manually by setting maxiter=1.
# -----------------------------
def validate_gradient_fd(md_in, index=2, delta=0.001):
    """
    Finite-difference check for dJ/dB at a given index (0-based python index),
    mirroring the MATLAB block. Requires md_in to already have outputdefinition
    and dependents set, and typically md_in.inversion.maxiter = 1.
    """
    md2 = md_in  # shallow reference; copy if you want isolation

    dJdB_ad = md2.results.TransientSolution[0].Gradient1[index]

    B1 = md2.materials.rheology_B[index].copy()
    B0 = B1 * (1.0 - delta)
    B2 = B1 * (1.0 + delta)
    deltaB = (B2 - B0)

    # requested outputs list from dependents
    out_list = [dep.name for dep in md2.autodiff.dependents]
    md2.transient.requested_outputs = out_list

    # turn off autodiff/inversion for pure forward evaluations
    md2.autodiff.isautodiff = False
    md2.inversion.iscontrol = False

    # forward at B0
    mdA = md2
    mdA.materials.rheology_B[index] = B0
    mdA = pyissm.model.execute.solve(mdA, 'tr')
    J0 = 0.0
    last = mdA.results.TransientSolution[-1]
    for nm in out_list:
        J0 += getattr(last, nm)

    # forward at B2
    mdB_ = md2
    mdB_.materials.rheology_B[index] = B2
    mdB_ = pyissm.model.execute.solve(mdB_, 'tr')
    J2 = 0.0
    last = mdB_.results.TransientSolution[-1]
    for nm in out_list:
        J2 += getattr(last, nm)

    dJdB_fd = (J2 - J0) / deltaB
    return dJdB_fd, dJdB_ad
