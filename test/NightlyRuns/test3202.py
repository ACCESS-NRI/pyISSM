# Test Name: SquareShelfTransientCalibrationNBVcodipack

import numpy as np
import pyissm

# -----------------------------
# Generate observations
# -----------------------------
md = pyissm.model.Model()
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 100000.)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelf.py')
md = pyissm.model.param.set_flow_equation(md, SSA='all')
md.cluster.np = 2

# -----------------------------
# Create real time series for B (element-based)
# -----------------------------
nv = md.mesh.numberofvertices
md.timestepping.interp_forcing = 0
md.timestepping.final_time = 2 * md.timestepping.time_step
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

# Forward solve to generate "observations"
md = pyissm.model.execute.solve(md, 'tr')

# -----------------------------
# Modify rheology, now constant
# -----------------------------
md.materials.rheology_B[:-1, :] = 1.8e8

# -----------------------------
# Set cost function (LogVel + Vy/Vx/Surface squares at each time)
# -----------------------------
md.outputdefinition.definitions = []
md.autodiff.dependents = []

count = 1
for i in range(0, len(md.results.TransientSolution.steps)):
    sol = md.results.TransientSolution[i]
    vx_obs = sol.Vx
    vy_obs = sol.Vy
    z_obs  = sol.Surface
    time   = sol.time
    weights = np.ones(nv)

    # 1) LogVel misfit
    od1 = pyissm.model.classes.cf.surfacelogvel()
    od1.name = f'LogVelMis{count}'
    od1.definitionstring = f'Outputdefinition{count}'
    od1.vxobs_string = 'VxObs'
    od1.vxobs = vx_obs
    od1.vyobs_string = 'VyObs'
    od1.vyobs = vy_obs
    od1.weights = weights
    od1.weights_string = 'WeightsSurfaceObservation'
    od1.datatime = time
    md.outputdefinition.definitions.append(od1)

    dep1 = pyissm.model.classes.dependent()
    dep1.name = f'Outputdefinition{count}'
    dep1.type = 'scalar'
    dep1.fos_reverse_index = 1
    dep1.nods = md.mesh.numberofvertices
    md.autodiff.dependents.append(dep1)
    count += 1

    # 2) Vy square misfit
    od2 = pyissm.model.classes.cf.surfacesquare()
    od2.name = f'VyMisfit{count}'
    od2.definitionstring = f'Outputdefinition{count}'
    od2.model_string = 'Vy'
    od2.observation_string = 'VyObs'
    od2.observation = vy_obs / md.constants.yts
    od2.weights = weights
    od2.weights_string = 'WeightsSurfaceObservation'
    od2.datatime = time
    md.outputdefinition.definitions.append(od2)

    dep2 = pyissm.model.classes.dependent()
    dep2.name = f'Outputdefinition{count}'
    dep2.type = 'scalar'
    dep2.fos_reverse_index = 1
    dep2.nods = md.mesh.numberofvertices
    md.autodiff.dependents.append(dep2)
    count += 1

    # 3) Vx square misfit (500*weights)
    od3 = pyissm.model.classes.cf.surfacesquare()
    od3.name = f'VxMisfit{count}'
    od3.definitionstring = f'Outputdefinition{count}'
    od3.model_string = 'Vx'
    od3.observation_string = 'VxObs'
    od3.observation = vx_obs / md.constants.yts
    od3.weights = 500.0 * weights
    od3.weights_string = 'WeightsSurfaceObservation'
    od3.datatime = time
    md.outputdefinition.definitions.append(od3)

    dep3 = pyissm.model.classes.dependent()
    dep3.name = f'Outputdefinition{count}'
    dep3.type = 'scalar'
    dep3.fos_reverse_index = 1
    dep3.nods = md.mesh.numberofvertices
    md.autodiff.dependents.append(dep3)
    count += 1

    # 4) DEM / surface misfit
    od4 = pyissm.model.classes.cf.surfacesquare()
    od4.name = f'DEMMisfit{count}'
    od4.definitionstring = f'Outputdefinition{count}'
    od4.model_string = 'Surface'
    od4.observation_string = 'SurfaceObservation'
    od4.observation = z_obs
    od4.weights = (1.0 / md.constants.yts) * weights
    od4.weights_string = 'WeightsSurfaceObservation'
    od4.datatime = time
    md.outputdefinition.definitions.append(od4)

    dep4 = pyissm.model.classes.dependent()
    dep4.name = f'Outputdefinition{count}'
    dep4.type = 'scalar'
    dep4.fos_reverse_index = 1
    dep4.nods = md.mesh.numberofvertices
    md.autodiff.dependents.append(dep4)
    count += 1

# -----------------------------
# Independent (control) on MaterialsRheologyBbar
# -----------------------------
min_params = md.materials.rheology_B.copy()
max_params = md.materials.rheology_B.copy()

# Replace these if cuffey is exposed elsewhere in your pyISSM.
min_params[:-1, :] = pyissm.tools.materials.cuffey(273)
max_params[:-1, :] = pyissm.tools.materials.cuffey(200)

indep = pyissm.model.classes.independent()
indep.name = 'MaterialsRheologyBbar'
indep.control_size = md.materials.rheology_B.shape[1]
indep.type = 'vertex'  # matches MATLAB comment (even if odd)
indep.min_parameters = min_params
indep.max_parameters = max_params
indep.control_scaling_factor = 1e8

md.autodiff.independents = [indep]

# -----------------------------
# Inversion + autodiff settings
# -----------------------------
md.inversion = pyissm.model.classes.inversion.adm1qn3(md.inversion)
md.inversion.iscontrol = 1
md.inversion.maxiter = 4
md.inversion.maxsteps = md.inversion.maxiter
md.inversion.dxmin = 1e-5

md.autodiff.isautodiff = True
md.autodiff.driver = 'fos_reverse'

# Go solve!
md = pyissm.model.execute.solve(md, 'tr')

# -----------------------------
# Fields and tolerances to track changes
# -----------------------------
field_names = ['Gradient', 'Misfit', 'Rheology']
field_tolerances = [2e-12, 1e-12, 1e-12]

ts1 = md.results.TransientSolution[0]  # MATLAB (1)
field_values = [
    ts1.Gradient1,
    ts1.J,
    ts1.MaterialsRheologyBbar
]
