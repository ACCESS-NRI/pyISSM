
#Test Name: PigStressMOLHO2d
import pyissm


md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Pig.exp', 20000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/PigShelves.exp', '../assets/Exp/PigIslands.exp')
md = pyissm.model.param.parameterize(md, '../assets/Par/Pig.py')
md = pyissm.model.param.set_flow_equation(md, MOLHO =  'all')
md.cluster.np = 3 
md.stressbalance.requested_outputs = ['default', 'VxSurface', 'VySurface', 'VxShear', 'VyShear', 'VxBase', 'VyBase']
md = pyissm.model.bc.set_molho_bc(md)
md = pyissm.model.execute.solve(md, 'Stressbalance')

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vel', 'Pressure', 'VxSurface', 'VySurface', 'VxShear', 'VyShear', 'VxBase', 'VyBase']
field_tolerances = [5e-13, 6e-13, 6e-13, 1e-13, 5e-13, 6e-13, 1e-13, 1e-13, 5e-13, 6e-13]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure,
                md.results.StressbalanceSolution.VxSurface,
                md.results.StressbalanceSolution.VySurface,
                md.results.StressbalanceSolution.VxShear,
                md.results.StressbalanceSolution.VyShear,
                md.results.StressbalanceSolution.VxBase,
                md.results.StressbalanceSolution.VyBase]
