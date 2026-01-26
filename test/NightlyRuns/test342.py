#Test Name: SquareSheetTherSteaPlume
import pyissm
import numpy as np

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 180000.)
md = pyissm.model.param.set_mask(md, None, None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetConstrained.py')
md.basalforcings = pyissm.model.classes.basalforcings.plume()
md.basalforcings.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices, ))
md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices, ))
md.basalforcings.plumex = 500000
md.basalforcings.plumey = 500000
md = md.extrude(3, 1.)
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3
md.timestepping.time_step = 0.
md.thermal.requested_outputs = ['default', 'BasalforcingsGeothermalflux']

# Execute model
md = pyissm.model.execute.solve(md, 'Thermal')

# Fields and tolerances to track changes
field_names = ['Temperature', 'BasalforcingsGroundediceMeltingRate', 'BasalforcingsGeothermalflux']
field_tolerances = [1e-13, 1e-8, 1e-13]
field_values = [md.results.ThermalSolution.Temperature,
                md.results.ThermalSolution.BasalforcingsGroundediceMeltingRate,
                md.results.ThermalSolution.BasalforcingsGeothermalflux]
