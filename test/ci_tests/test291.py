#Test Name: SquareShelfStressFSP4z
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelf.py')
md = md.extrude(2, 1.)
md = pyissm.model.param.set_flow_equation(md, FS = 'all')
md.flowequation.fe_FS = 'OneLayerP4z'
md.cluster.np = 3

# Execute model
md = pyissm.model.execute.solve(md, 'Stressbalance')

# Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure']
field_tolerances = [5e-5, 5e-5, 1e-6, 9e-5, 1e-7]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vz,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure]
