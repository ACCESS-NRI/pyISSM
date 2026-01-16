#Test Name: SquareShelfTherStea
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 180000)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelf.py')
md = md.extrude(3, 1.)
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.timestepping.time_step = 0
md.cluster.np = 3

# Execute model
md = pyissm.model.execute.solve(md, 'Thermal')

#  Fields and tolerances to track changes
field_names = ['Temperature', 'BasalforcingsGroundediceMeltingRate']
field_tolerances = [1e-13, 6e-6]
field_values = [md.results.ThermalSolution.Temperature, md.results.ThermalSolution.BasalforcingsGroundediceMeltingRate]
