#Test Name: PigTherSteaSUPG
import pyissm

md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Pig.exp', 30000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/PigShelves.exp', '../assets/Exp/PigIslands.exp')
md = pyissm.model.param.parameterize(md, '../assets/Par/Pig.py')
md = md.extrude(3, 1.)
md = pyissm.model.setflowequation(md, HO = 'all')
md.thermal.stabilization = 2
md.cluster.np = 3 
md.timestepping.time_step = 0
md.thermal.penalty_threshold = 40
md = pyissm.model.execute.solve(md, 'Thermal')
#Fields and tolerances to track changes
field_names = ['Temperature', 'BasalforcingsGroundediceMeltingRate']
field_tolerances = [1e-11, 1e-11]
field_values = [md.results.ThermalSolution.Temperature,
                md.results.ThermalSolution.BasalforcingsGroundediceMeltingRate]
