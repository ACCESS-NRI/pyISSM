#Test Name: PigSteaHO
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Pig.exp', 30000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/PigShelves.exp', '../assets/Exp/PigIslands.exp')
md = pyissm.model.param.parameterize(md, '../assets/Par/Pig.py')
md = md.extrude(3, 1.)
md = pyissm.model.param.set_flow_equation(md, HO = 'all')
md.cluster.np = 3
md.timestepping.time_step = 0.
md.thermal.penalty_threshold = 7

# Execute model
md = pyissm.model.execute.solve(md, 'Steadystate')

# Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure', 'Temperature', 'BasalforcingsGroundediceMeltingRate']
field_tolerances = [3e-09, 3e-09, 5e-08, 5e-08, 1e-09, 8e-09, 1e-06]
field_values = [md.results.SteadystateSolution.Vx,
                md.results.SteadystateSolution.Vy,
                md.results.SteadystateSolution.Vz,
                md.results.SteadystateSolution.Vel,
                md.results.SteadystateSolution.Pressure,
                md.results.SteadystateSolution.Temperature,
                md.results.SteadystateSolution.BasalforcingsGroundediceMeltingRate]
