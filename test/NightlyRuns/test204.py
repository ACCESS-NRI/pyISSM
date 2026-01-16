#Test Name: SquareShelfStressFS
import pyissm
import copy

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 180000)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelf.py')
md = md.extrude(3, 2.)
md = pyissm.model.param.set_flow_equation(md, FS = 'all')
md.cluster.np = 1
md.stressbalance.shelf_dampening = 1
md.timestepping.time_step = 0

# Execute model
md1 = copy.deepcopy(md)
md1 = pyissm.model.execute.solve(md1, 'Stressbalance')

# Execute second model (shelf_dampening = 0)
md.stressbalance.shelf_dampening = 0
md = pyissm.model.execute.solve(md, 'Stressbalance')

#  Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure','Vx_damp','Vy_damp','Vz_damp','Vel_damp','Pressure_damp']
field_tolerances = [1e-08, 1e-08, 8e-06, 1e-08, 1e-08, 1e-08, 1e-08, 2e-07, 1e-08, 1e-08]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vz,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure,
                md1.results.StressbalanceSolution.Vx,
                md1.results.StressbalanceSolution.Vy,
                md1.results.StressbalanceSolution.Vz,
                md1.results.StressbalanceSolution.Vel,
                md1.results.StressbalanceSolution.Pressure]