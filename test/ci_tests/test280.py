#Test Name: SquareShelfStressSSA2dHigherOrder
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelf.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3

# Execute model and define fields and tolerances to track changes
field_names = []
field_tolerances = []
field_values = []
for i in ['P1bubble', 'P1bubblecondensed', 'P2']:
    md.flowequation.fe_SSA = i
    md = pyissm.model.execute.solve(md, 'Stressbalance')
    field_names = field_names + ['Vx' + i, 'Vy' + i, 'Vel' + i, 'Pressure' + i]
    field_tolerances = field_tolerances + [1e-12, 1e-13, 1e-13, 1e-13]
    field_values = field_values + [md.results.StressbalanceSolution.Vx,
                                   md.results.StressbalanceSolution.Vy,
                                   md.results.StressbalanceSolution.Vel,
                                   md.results.StressbalanceSolution.Pressure]