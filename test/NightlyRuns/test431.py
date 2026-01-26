#Test Name: SquareSheetShelfSteaEnthalpySSA3d
import pyissm
import numpy as np

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/SquareShelf.exp', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetShelf.py')
md = md.extrude(3, 2.)
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3
md.timestepping.time_step = 0.
md.thermal.isenthalpy = 1
md.thermal.isdynamicbasalspc = 1
md.initialization.waterfraction = np.zeros((md.mesh.numberofvertices))
md.initialization.watercolumn = np.zeros((md.mesh.numberofvertices))

# Execute model
md = pyissm.model.execute.solve(md, 'Steadystate')

# Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure', 'Temperature', 'Waterfraction', 'Enthalpy']
field_tolerances = [8e-4, 5e-4, 5e-4, 5e-4, 1e-13, 1e-4, 5e-4, 5e-4]
field_values = [md.results.SteadystateSolution.Vx,
                md.results.SteadystateSolution.Vy,
                md.results.SteadystateSolution.Vz,
                md.results.SteadystateSolution.Vel,
                md.results.SteadystateSolution.Pressure,
                md.results.SteadystateSolution.Temperature,
                md.results.SteadystateSolution.Waterfraction,
                md.results.SteadystateSolution.Enthalpy]
