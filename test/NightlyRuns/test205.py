#Test Name: SquareShelfStressMHOPenalties
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelf.py')
md = md.extrude(3, 2.)
md = pyissm.model.param.set_flow_equation(md,
                                          HO = '../assets/Exp/SquareHalfRight.exp',
                                          fill = 'SSA',
                                          coupling = 'penalties')
md.cluster.np = 3
md.settings.solver_residue_threshold = 1.e-4

# Execute model
md = pyissm.model.execute.solve(md, 'Stressbalance')


# Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure']
field_tolerances = [2e-05, 2e-05, 1e-05, 1e-05, 1e-05]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vz,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure]
