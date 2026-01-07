#Test Name: SquareShelfStressSSA2d
import pyissm

md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelf.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3

# TEMPORARY FIX UNTIL WAITONLOCK IS FIXED
md.settings.waitonlock = 0
md.miscellaneous.name = 'SquareShelfStressSSA2d'
md = pyissm.model.execute.solve(md, 'Stressbalance')

# ## TEMPORARY FIX UNTIL WAITONLOCK IS FIXED
# # Once fixed, remove runtime_name=False from above and below and load_only=True from below. Remove waitonlock = 0 above. This means python will then wait for the results to be written to disk before loading them.
# import time
# print('Waiting for results to be written to disk...')
# time.sleep(5)  #wait for results, in absence of waitonlock
# print('Loading results from disk...')
# md = pyissm.model.execute.solve(md, 'Stressbalance', runtime_name=False, load_only=True)
# print(md.results.StressbalanceSolution)

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vel', 'Pressure']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure]
