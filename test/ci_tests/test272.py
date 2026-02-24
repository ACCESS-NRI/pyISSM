#Test Name: SquareShelfCMDSSA2dDamage
import numpy as np
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000)
md = pyissm.model.param.set_mask(md, 'all', None)
md.materials = pyissm.model.classes.materials.damageice()
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelf.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3
md.damage.isdamage = 1
md.damage.D = 0.5 * np.ones(md.mesh.numberofvertices)
md.damage.spcdamage = np.nan * np.ones(md.mesh.numberofvertices)

# control parameters
md.inversion.iscontrol = 1
md.inversion.control_parameters = ['DamageDbar']
md.inversion.iscontrol = 1
md.inversion.control_parameters = ['DamageDbar']
md.inversion.min_parameters = 10**-13 * np.ones((md.mesh.numberofvertices, len(md.inversion.control_parameters)))
md.inversion.max_parameters = np.ones((md.mesh.numberofvertices, len(md.inversion.control_parameters)))
md.inversion.nsteps = 2
md.inversion.cost_functions = [101]
md.inversion.cost_functions_coefficients = np.ones((md.mesh.numberofvertices, len(md.inversion.cost_functions)))
md.inversion.gradient_scaling = 0.9 * np.ones((md.inversion.nsteps, len(md.inversion.control_parameters)))
md.inversion.maxiter_per_step = 2. * np.ones((md.inversion.nsteps))
md.inversion.step_threshold = 0.99 * np.ones((md.inversion.nsteps))
md.inversion.vx_obs = md.initialization.vx
md.inversion.vy_obs = md.initialization.vy

# Execute model
md = pyissm.model.execute.solve(md, 'Stressbalance')

# Fields and tolerances to track changes
field_names = ['Gradient', 'Misfits', 'DamageDbar', 'Pressure', 'Vel', 'Vx', 'Vy']
field_tolerances = [1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12]
field_values = [md.results.StressbalanceSolution.Gradient1,
                md.results.StressbalanceSolution.J,
                md.results.StressbalanceSolution.DamageDbar,
                md.results.StressbalanceSolution.Pressure,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy]
