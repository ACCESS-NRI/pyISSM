# Test Name: SquareShelfAdolcStaticControls
import numpy as np
import pyissm

# --- mesh + mask + params ---
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../Exp/Square.exp', 100000.)
md = pyissm.model.param.set_mask(md, '../Exp/SquareShelf.exp', '')   # same intent as setmask(exp,'')

md = pyissm.model.param.parameterize(md, '../Par/SquareSheetShelf.par')

# --- initialization ---
md.initialization.vx[:] = 1.0
md.initialization.vy[:] = 1.0

# thickness = 500 - x/10000
md.geometry.thickness[:] = 500.0 - md.mesh.x / 10000.0

# bed = -100 - x/1000
md.geometry.bed[:] = -100.0 - md.mesh.x / 1000.0

# base = -H * rho_i / rho_w
md.geometry.base[:] = -md.geometry.thickness * md.materials.rho_ice / md.materials.rho_water

# ocean_levelset = H + (rho_w/rho_i)*bed
md.mask.ocean_levelset[:] = md.geometry.thickness + (md.materials.rho_water / md.materials.rho_ice) * md.geometry.bed

# pos = find(ocean_levelset>=0); base(pos)=bed(pos)
pos = np.where(md.mask.ocean_levelset >= 0)[0]
md.geometry.base[pos] = md.geometry.bed[pos]

# surface = base + thickness
md.geometry.surface[:] = md.geometry.base + md.geometry.thickness

# --- flow equation ---
md = pyissm.model.param.set_flow_equation(md, SSA='all')

# --- control parameters / inversion / autodiff ---
md.inversion = pyissm.model.classes.adm1qn3inversion(md.inversion)
md.inversion.iscontrol = 1

md.autodiff.isautodiff = True
md.autodiff.driver = 'fos_reverse'

# friction coefficient init
md.friction.coefficient[:md.mesh.numberofvertices] = 50.0

# independents: FrictionCoefficient control (vertex)
indep = pyissm.model.classes.independent()
indep.name = 'FrictionCoefficient'
indep.type = 'vertex'
indep.control_size = 1
indep.min_parameters = np.ones((md.mesh.numberofvertices,))
indep.max_parameters = 500.0 * np.ones((md.mesh.numberofvertices,))
indep.control_scaling_factor = 1

md.autodiff.independents = [indep]

# --- outputdefinition: cfsurfacesquare ---
# observation = vy/yts
vy_obs = md.initialization.vy / md.constants.yts
weights = np.ones((md.mesh.numberofvertices,))

od = pyissm.model.classes.cfsurfacesquare()
od.name = 'VyMisfit1'
od.definitionstring = 'Outputdefinition1'
od.model_string = 'Vy'
od.observation_string = 'VyObs'
od.observation = vy_obs
od.weights = weights
od.weights_string = 'WeightsSurfaceObservation'
od.datatime = 0.75

# put into definitions{1}
md.outputdefinition.definitions = [od]

# --- timestepping / transient flags ---
md.timestepping.interp_forcing = 0
md.timestepping.time_step = 0.5
md.timestepping.final_time = 1.5

md.transient.ismasstransport = 1
md.transient.isstressbalance = 1
md.transient.isgroundingline = 1
md.transient.ismovingfront = 0
md.transient.isthermal = 0

# basal forcings
md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices,))
md.basalforcings.floatingice_melting_rate = 25.0 * np.zeros((md.mesh.numberofvertices,))

# --- dependents: Outputdefinition1 scalar with reverse index 1 ---
dep = pyissm.model.classes.dependent()
dep.name = 'Outputdefinition1'
dep.type = 'scalar'
dep.fos_reverse_index = 1

md.autodiff.dependents = [dep]

# inversion iterations
md.inversion.maxiter = 2

# cluster
md.cluster = pyissm.cluster.generic('name', pyissm.tools.oshostname(), 'np', 3)

# --- solve ---
md = pyissm.execute.solve(md, 'transient')

# --- Fields and tolerances to track changes ---
field_names = [
    'Gradient', 'Misfits', 'FrictionCoefficient',
    'Pressure1', 'Vel1', 'Vx1', 'Vy1',
    'Pressure2', 'Vel2', 'Vx2', 'Vy2'
]

field_tolerances = [1e-11] * len(field_names)

# MATLAB used TransientSolution(1) and (3).
# In python, results are typically 0-indexed, so (1)->0 and (3)->2
ts1 = md.results.TransientSolution[0]
ts3 = md.results.TransientSolution[2]

field_values = [
    ts1.Gradient1,
    ts1.J,
    ts1.FrictionCoefficient,
    ts1.Pressure,
    ts1.Vel,
    ts1.Vx,
    ts1.Vy,
    ts3.Pressure,
    ts3.Vel,
    ts3.Vx,
    ts3.Vy
]
