#Test Name: SquareShelfConstrainedTherTranAdolc
import pyissm


md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 180000.)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
md = md.extrude(3, 1.)
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 1
md.transient.isstressbalance = False
md.transient.ismasstransport = False
md.transient.issmb = True
md.transient.isthermal = True
md.transient.isgroundingline = False
md.autodiff.isautodiff = True
md.toolkits.DefaultAnalysis = pyissm.tools.config.issm_gsl_solver()
md = pyissm.model.execute.solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Temperature', 'BasalforcingsGroundediceMeltingRate']
field_tolerances = [1e-13, 1e-13]
field_values = [md.results.TransientSolution[0].Temperature,
                md.results.TransientSolution[0].BasalforcingsGroundediceMeltingRate]
