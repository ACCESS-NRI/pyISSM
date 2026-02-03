#Test Name: SquareShelfConstrainedMasstransp3dAdolcMumps
import pyissm


md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000.)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md = md.extrude(5, 3.)
md.cluster.np = 3

md.autodiff.isautodiff = True
md.toolkits.DefaultAnalysis = pyissm.tools.config.issm_mumps_solver()
md = pyissm.model.execute.solve(md, 'Masstransport')

#Fields and tolerances to track changes
field_names = ['Thickness']
field_tolerances = [1e-13]
field_values = [md.results.MasstransportSolution.Thickness]
