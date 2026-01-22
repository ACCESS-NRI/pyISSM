#Test Name: SquareShelfConstrainedEnthalpyStea
import pyissm

md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 180000)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
md = md.extrude(3, 1.)
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.timestepping.time_step = 0
md.initialization.waterfraction = np.zeros(md.mesh.numberofvertices)
md.initialization.watercolumn = np.zeros(md.mesh.numberofvertices)
md.thermal.isenthalpy = 1
md.thermal.isdynamicbasalspc = 1

md.cluster.np = 3
md = pyissm.model.execute.solve(md, 'Thermal')

#Fields and tolerances to track changes
field_names = ['Enthalpy', 'Waterfraction', 'Temperature']
field_tolerances = [1e-13, 2e-10, 1e-13]
field_values = [md.results.ThermalSolution.Enthalpy,
                md.results.ThermalSolution.Waterfraction,
                md.results.ThermalSolution.Temperature]
