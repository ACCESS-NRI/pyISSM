#Test Name: SquareSheetConstrainedStressHO
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 180000.)
md = pyissm.model.param.set_mask(md, None, None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetConstrained.py')
md = md.extrude(3, 1.)
md = pyissm.model.param.set_flow_equation(md, HO = 'all')
md.cluster.np = 3

# Execute model
md = pyissm.model.execute.solve(md, 'Stressbalance')

# Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure']
field_tolerances = [1e-09, 1e-09, 2e-10, 2e-10, 1e-10]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vz,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure]
