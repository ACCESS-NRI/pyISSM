#Test Name: SquareSheetConstrainedTherTranNyeH2O
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 180000.)
md = pyissm.model.param.set_mask(md, None, None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetConstrained.py')
md = md.extrude(3, 1.)
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3

md.materials.rheology_law = 'NyeH2O'
md.materials.rheology_B = pyissm.tools.materials.nye(md.initialization.temperature, 2)

# Transient options
md.transient.isstressbalance = 0
md.transient.ismasstransport = 0
md.transient.issmb = 1
md.transient.isthermal = 1
md.transient.isgroundingline = 0

# Execute model
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Temperature1', 'BasalforcingsGroundediceMeltingRate1',
               'Temperature3', 'BasalforcingsGroundediceMeltingRate3']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md.results.TransientSolution[0].Temperature,
                md.results.TransientSolution[0].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[2].Temperature,
                md.results.TransientSolution[2].BasalforcingsGroundediceMeltingRate]
