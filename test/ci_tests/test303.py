#Test Name: SquareSheetConstrainedStressSIASSA2d
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000.)
md = pyissm.model.param.set_mask(md, None, None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetConstrained.py')
md = pyissm.model.param.set_flow_equation(md, SIA = '../assets/Exp/SquareHalfRight.exp', fill = 'SSA')
md.cluster.np = 3

# Execute model
md = pyissm.model.execute.solve(md, 'Stressbalance')

# Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vel', 'Pressure']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure]
