#Test Name: SquareSheetConstrainedMasstransp2dDG
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000.)
md = pyissm.model.mesh.mesh_convert(md)
md = pyissm.model.param.set_mask(md, None, None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetConstrained.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.masstransport.stabilization = 3
md.masstransport.spcthickness = md.geometry.thickness
md.cluster.np = 3

# Excute model
md = pyissm.model.execute.solve(md, 'Masstransport')

# Fields and tolerances to track changes
field_names = ['Thickness']
field_tolerances = [1e-13]
field_values = [md.results.MasstransportSolution.Thickness]
