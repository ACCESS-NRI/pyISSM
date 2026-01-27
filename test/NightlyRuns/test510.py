#Test Name: PigSteaFS
import pyissm


md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Pig.exp', 20000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/PigShelves.exp', '../assets/Exp/PigIslands.exp')
md = pyissm.model.param.parameterize(md, '../assets/Par/Pig.py')
md = md.extrude(2, 1.)
md = pyissm.model.param.set_flow_equation(md, FS = 'all')
md.cluster.np = 3
md.timestepping.time_step = 0.
md = pyissm.model.execute.solve(md, 'Steadystate')

# Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure', 'Temperature', 'BasalforcingsGroundediceMeltingRate']
field_tolerances = [1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-06]
field_values = [md.results.SteadystateSolution.Vx,
                md.results.SteadystateSolution.Vy,
                md.results.SteadystateSolution.Vz,
                md.results.SteadystateSolution.Vel,
                md.results.SteadystateSolution.Pressure,
                md.results.SteadystateSolution.Temperature,
                md.results.SteadystateSolution.BasalforcingsGroundediceMeltingRate]
