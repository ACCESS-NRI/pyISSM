#Test Name: SquareSheetShelfStressSSAFS3dTiling
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 180000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/SquareShelf.exp', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetShelf.py')
md = md.extrude(5, 1.)
md = pyissm.model.param.set_flow_equation(md, FS = '../assets/Exp/SquareHalfRight.exp', fill = 'SSA')
md.cluster.np = 3
md.stressbalance.reltol = 0.4

# Execute model
md = pyissm.model.execute.solve(md, 'Stressbalance')

# Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure']
field_tolerances = [4e-07, 4e-07, 2e-06, 4e-07, 5e-07]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vz,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure]
