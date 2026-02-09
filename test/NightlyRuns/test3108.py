#Test Name: SquareShelfConstrainedTherSteaAdolcMumps
import pyissm


md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 180000.)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
md = md.extrude(3, 1.)
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.timestepping.time_step = 0
md.cluster.np = 3
md.autodiff.isautodiff = True
md.toolkits.DefaultAnalysis = pyissm.tools.config.issm_mumps_solver()
md = pyissm.model.execute.solve(md, 'Thermal')

#Fields and tolerances to track changes
field_names = ['Temperature', 'BasalforcingsGroundediceMeltingRate']
field_tolerances = [1e-13, 1e-5]
field_values = [md.results.ThermalSolution.Temperature,
                md.results.ThermalSolution.BasalforcingsGroundediceMeltingRate]
