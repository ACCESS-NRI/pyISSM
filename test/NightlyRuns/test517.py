#Test Name: PigCMBFSm1qn3
import pyissm
import numpy as np

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Pig.exp', 11000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/PigShelves.exp', '../assets/Exp/PigIslands.exp')
md = pyissm.model.param.parameterize(md, '../assets/Par/Pig.py')

# impose hydrostatic equilibrium (required by Stokes)
md.geometry.base = -md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness
md.geometry.surface = md.geometry.base + md.geometry.thickness
md = md.extrude(3, 1.)
md = pyissm.model.param.set_flow_equation(md, FS = 'all')
md = md.extract(md.mask.ocean_levelset < 0.)

# control parameters
md.inversion.iscontrol = 1
md.inversion.control_parameters = ['MaterialsRheologyBbar']
md.inversion.min_parameters = 10.**6 * np.ones((md.mesh.numberofvertices, len(md.inversion.control_parameters)))
md.inversion.max_parameters = 2. * 10**9 * np.ones((md.mesh.numberofvertices, len(md.inversion.control_parameters)))
md.inversion.nsteps = 2
md.inversion.cost_functions = [101, 502]
md.inversion.cost_functions_coefficients = np.ones((md.mesh.numberofvertices, len(md.inversion.cost_functions)))
md.inversion.gradient_scaling = 10.**8 * np.ones((md.inversion.nsteps, len(md.inversion.control_parameters)))
md.inversion.maxiter_per_step = 2. * np.ones((md.inversion.nsteps))
md.inversion.step_threshold = 0.99 * np.ones((md.inversion.nsteps))
md.inversion.vx_obs = md.initialization.vx
md.inversion.vy_obs = md.initialization.vy

md.inversion = pyissm.model.classes.inversion.m1qn3(md.inversion)

# Execute model
md = pyissm.model.execute.solve(md, 'Stressbalance')

# Fields and tolerances to track changes
field_names = ['Gradient', 'Misfits', 'MaterialsRheologyB', 'Pressure', 'Vel', 'Vx', 'Vy']
field_tolerances = [6e-11, 5e-11, 5e-10, 1e-09, 2e-11, 5e-11, 2e-11]
field_values = [md.results.StressbalanceSolution.Gradient1,
                md.results.StressbalanceSolution.J,
                md.results.StressbalanceSolution.MaterialsRheologyBbar,
                md.results.StressbalanceSolution.Pressure,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy]
