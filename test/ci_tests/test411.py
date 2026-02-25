#Test Name: SquareSheetShelfSteaHO
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/SquareShelf.exp', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetShelf.py')
md = md.extrude(3, 2.)
md = pyissm.model.param.set_flow_equation(md, HO = 'all')
md.cluster.np = 3
md.timestepping.time_step = 0

# Execute model
md = pyissm.model.execute.solve(md, 'Steadystate')

# Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure', 'Temperature', 'BasalforcingsGroundediceMeltingRate']
field_tolerances = [2e-09, 1e-09, 1e-08, 1e-09, 1e-09, 1e-08, 1e-05]
field_values = [md.results.SteadystateSolution.Vx,
                md.results.SteadystateSolution.Vy,
                md.results.SteadystateSolution.Vz,
                md.results.SteadystateSolution.Vel,
                md.results.SteadystateSolution.Pressure,
                md.results.SteadystateSolution.Temperature,
                md.results.SteadystateSolution.BasalforcingsGroundediceMeltingRate]
