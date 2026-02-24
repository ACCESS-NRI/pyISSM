#Test Name: SquareShelfTranSSA3d
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 180000)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelf.py')
md = md.extrude(3, 1.)
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3

# Execute model
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vz1', 'Vel1', 'Pressure1', 'Bed1', 'Surface1', 'Thickness1', 'Temperature1', 'BasalforcingsGroundediceMeltingRate1',
               'Vx2', 'Vy2', 'Vz2', 'Vel2', 'Pressure2', 'Bed2', 'Surface2', 'Thickness2', 'Temperature2', 'BasalforcingsGroundediceMeltingRate2',
               'Vx3', 'Vy3', 'Vz3', 'Vel3', 'Pressure3', 'Bed3', 'Surface3', 'Thickness3', 'Temperature3', 'BasalforcingsGroundediceMeltingRate3']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-8,
                    1e-13, 1e-13, 3e-9, 1e-13, 3.8e-11, 3.85e-11, 3.85e-11, 3.85e-11, 1e-13, 3e-8,
                    5e-11, 6e-12, 1e-08, 6e-12, 3.9e-11, 1e-10, 1e-10, 1e-10, 9e-12, 5e-8]

field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vz,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].Temperature,
                md.results.TransientSolution[0].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Vz,
                md.results.TransientSolution[1].Vel,
                md.results.TransientSolution[1].Pressure,
                md.results.TransientSolution[1].Base,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].Temperature,
                md.results.TransientSolution[1].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Vz,
                md.results.TransientSolution[2].Vel,
                md.results.TransientSolution[2].Pressure,
                md.results.TransientSolution[2].Base,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].Temperature,
                md.results.TransientSolution[2].BasalforcingsGroundediceMeltingRate]
