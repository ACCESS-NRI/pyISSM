#Test Name: SquareShelfConstrainedStressSSA3d
import pyissm


md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 180000)
md = pyissm.model.param.set_mask(md, 'all', '')
md = pyissm.model.param.parameterize(md, '../assets/ar/SquareShelfConstrained.py')
md.extrude(3, 2.)
md = pyissm.model.param.set_flow_equation(md, 'SSA', 'all')
md.cluster.name = 3
md = pyissm.model.execute.solve(md, 'Stressbalance')
#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vz,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure]
