#Test Name: PigBalVel1
import pyissm

md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Pig.exp', 20000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/PigShelves.exp', '../assets/Exp/PigIslands.exp')
md = pyissm.model.param.parameterize(md, '../assets/Par/Pig.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3 
md = pyissm.model.execute.solve(md, 'balancevelocity')

# Fields and tolerances to track changes
field_names = ['DrivingStressX', 'DrivingStressY', 'Vel']
field_tolerances = [1e-13, 1e-13, 1e-13]
field_values = [md.results.BalancevelocitySolution.DrivingStressX,
                md.results.BalancevelocitySolution.DrivingStressY,
                md.results.BalancevelocitySolution.Vel]
