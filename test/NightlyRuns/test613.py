#Test Name: 79NorthCMBalThicVxVy
import pyissm
import numpy as np
import copy

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/79North.exp', 10000.)
md = pyissm.model.mesh.mesh_convert(md)
md = pyissm.model.param.set_mask(md, '../assets/Exp/79NorthShelf.exp', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/79North.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3

# Ice sheet only
md = md.extract(md.mask.ocean_levelset > 0.)
pos = np.nonzero(md.mesh.vertexonboundary)
md.balancethickness.spcthickness[pos] = md.geometry.thickness[pos]

# control parameters
md.inversion.thickness_obs = copy.deepcopy(md.geometry.thickness)
md.inversion.iscontrol = 1
md.inversion.nsteps = 2
md.inversion.control_parameters = ['Vx', 'Vy']
md.balancethickness.stabilization = 1
md.inversion.gradient_scaling = np.vstack((10. / md.constants.yts * np.ones((md.inversion.nsteps)), 10. / md.constants.yts * np.ones((md.inversion.nsteps)))).T
md.inversion.min_parameters = np.vstack((-2000. * np.ones((md.mesh.numberofvertices)), -2000. * np.ones((md.mesh.numberofvertices)))).T
md.inversion.max_parameters = np.vstack((+ 2000. * np.ones((md.mesh.numberofvertices)), +2000. * np.ones((md.mesh.numberofvertices)))).T
md.inversion.cost_functions = [201]
md.inversion.cost_functions_coefficients = np.ones((md.mesh.numberofvertices, len(md.inversion.cost_functions)))
md.inversion.maxiter_per_step = 4 * np.ones((md.inversion.nsteps))
md.inversion.step_threshold = 0.99 * np.ones((md.inversion.nsteps))

md.verbose.control = 1

# Execute model
md = pyissm.model.execute.solve(md, 'Balancethickness')

# Fields and tolerances to track changes
field_names = ['Gradient1', 'Gradient2', 'Misfits', 'Vx', 'Vy', 'Thickness']
field_tolerances = [1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12]
field_values = [md.results.BalancethicknessSolution.Gradient1,
                md.results.BalancethicknessSolution.Gradient2,
                md.results.BalancethicknessSolution.J,
                md.results.BalancethicknessSolution.Vx,
                md.results.BalancethicknessSolution.Vy,
                md.results.BalancethicknessSolution.Thickness]
