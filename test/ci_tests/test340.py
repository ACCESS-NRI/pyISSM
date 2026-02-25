#Test Name: SquareSheetConstrainedSmbGradientsEla3d
import pyissm
import numpy as np

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 200000.)
md = pyissm.model.param.set_mask(md, None, None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetConstrained.py')
md = md.extrude(3, 1.)
md = pyissm.model.param.set_flow_equation(md, HO = 'all')
md.cluster.np = 3

# control parameters
md.inversion = pyissm.model.classes.inversion.tao()
md.inversion.iscontrol = 1
md.inversion.control_parameters = ['FrictionCoefficient']
md.inversion.min_parameters = 1. * np.ones((md.mesh.numberofvertices, len(md.inversion.control_parameters)))
md.inversion.max_parameters = 200. * np.ones((md.mesh.numberofvertices, len(md.inversion.control_parameters)))
md.inversion.maxsteps = 2
md.inversion.maxiter = 6
md.inversion.cost_functions = [102, 501]
md.inversion.cost_functions_coefficients = np.ones((md.mesh.numberofvertices, 2))
md.inversion.cost_functions_coefficients[:, 1] = 2. * 1e-7
md.inversion.vx_obs = md.initialization.vx
md.inversion.vy_obs = md.initialization.vy

# Execute model
md = pyissm.model.execute.solve(md, 'Stressbalance')

# Fields and tolerances to track changes
field_names = ['Gradient', 'Misfits', 'FrictionCoefficient', 'Pressure', 'Vel', 'Vx', 'Vy']
field_tolerances = [3e-08, 1e-07, 5e-10, 1e-10, 1e-09, 1e-09, 1e-09]
field_values = [md.results.StressbalanceSolution.Gradient1,
                md.results.StressbalanceSolution.J,
                md.results.StressbalanceSolution.FrictionCoefficient,
                md.results.StressbalanceSolution.Pressure,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy]
