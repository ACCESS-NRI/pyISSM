#Test Name: PigTherTranSUPG
import pyissm


md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Pig.exp', 30000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/PigShelves.exp', '../assets/Exp/PigIslands.exp')
md = pyissm.model.param.parameterize(md, '../assets/Par/Pig.py')
md = md.extrude(3, 1.)
md = pyissm.model.setflowequation(md, HO = 'all')
md.thermal.stabilization = 2
md.cluster.np = 3 
md.transient.isstressbalance = False
md.transient.ismasstransport = False
md.transient.issmb = True
md.transient.isthermal = True
md.transient.isgroundingline = False
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Temperature1', 'BasalforcingsGroundediceMeltingRate1',
               'Temperature2', 'BasalforcingsGroundediceMeltingRate2']
field_tolerances = [1e-13, 1e-8, 1e-13, 5e-8]
field_values = [md.results.TransientSolution[0].Temperature,
                md.results.TransientSolution[0].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[1].Temperature,
                md.results.TransientSolution[1].BasalforcingsGroundediceMeltingRate]
