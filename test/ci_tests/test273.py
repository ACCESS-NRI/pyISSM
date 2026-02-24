#Test Name: SquareShelfStressSSA2dDamageUpdate
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
md.damage.D = 0. * np.ones(md.mesh.numberofvertices)
md.damage.spcdamage = np.nan * np.ones(md.mesh.numberofvertices)
md.stressbalance.requested_outputs = ['default', 'NewDamage']
md.damage.stress_threshold = 1.3e5
md.damage.kappa = 2.8

# Execute model
md = pyissm.model.execute.solve(md, 'Stressbalance')

# Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vel', 'Pressure', 'NewDamage']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure,
                md.results.StressbalanceSolution.NewDamage]
