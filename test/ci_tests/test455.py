#Test Name: SquareSheetShelfStressHOHigherOrder
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 180000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/SquareShelf.exp', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetShelf.py')
md = md.extrude(5, 1.)
md = pyissm.model.param.set_flow_equation(md, HO = 'all')
md.cluster.np = 3

# Execute model
field_names = []
field_tolerances = []
field_values = []
for i in ['P1bubble', 'P1bubblecondensed', 'P1xP2', 'P2xP1', 'P2', 'P1xP3', 'P2xP4']:
    md.flowequation.fe_HO = i
    md = pyissm.model.execute.solve(md, 'Stressbalance')
    field_names = field_names + ['Vx' + i, 'Vy' + i, 'Vz' + i, 'Vel' + i, 'Pressure' + i]
    field_tolerances = field_tolerances + [1e-07, 6e-08, 6e-08, 6e-08, 3e-13]
    field_values = field_values + [md.results.StressbalanceSolution.Vx,
                                   md.results.StressbalanceSolution.Vy,
                                   md.results.StressbalanceSolution.Vz,
                                   md.results.StressbalanceSolution.Vel,
                                   md.results.StressbalanceSolution.Pressure]
