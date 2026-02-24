#Test Name: SquareShelfStressSSAFS3dTiling
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 120000)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelf2.py')
md = md.extrude(2, 1.)
md = pyissm.model.param.set_flow_equation(md,
                                          FS = '../assets/Exp/SquareHalfRight.exp',
                                          fill = 'SSA')
md.cluster.np = 3

# Execute model
md = pyissm.model.execute.solve(md, 'Stressbalance')

# Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel']
field_tolerances = [1e-09, 1e-09, 5e-06, 1e-09]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vz,
                md.results.StressbalanceSolution.Vel
]
