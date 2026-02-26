#Test Name: SquareShelfConstrainedRestartTranHO3d
import pyissm
import copy

md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000.)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
md.cluster.np = 1
md.transient.requested_outputs = ['IceVolume', 'TotalSmb']
md = md.extrude(3, 1.)
md = pyissm.model.param.set_flow_equation(md, HO = 'all')

md.settings.checkpoint_frequency = 5

# time steps and resolution
md.timestepping.final_time = 8

md = pyissm.model.execute.solve(md, 'Transient')
md2 = copy.deepcopy(md)
md = pyissm.model.execute.solve(md, 'Transient', 'restart', 1)

#Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'TotalSmb1', 'Bed1', 'Surface1', 'Thickness1', 'Volume1', 'Temperature1', 'Pressure1', 'Vx2', 'Vy2', 'Vel2', 'TotalSmb2', 'Bed2', 'Surface2', 'Thickness2', 'Volume2', 'Temperature2', 'Pressure2', 'Vx3', 'Vy3', 'Vel3', 'TotalSmb3', 'Bed3', 'Surface3', 'Thickness3', 'Volume3', 'Temperature3', 'Pressure3']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md2.results.TransientSolution[5].Vx - md.results.TransientSolution[5].Vx,
                md2.results.TransientSolution[5].Vy - md.results.TransientSolution[5].Vy,
                md2.results.TransientSolution[5].Vel - md.results.TransientSolution[5].Vel,
                md2.results.TransientSolution[5].TotalSmb - md.results.TransientSolution[5].TotalSmb,
                md2.results.TransientSolution[5].Base - md.results.TransientSolution[5].Base,
                md2.results.TransientSolution[5].Surface - md.results.TransientSolution[5].Surface,
                md2.results.TransientSolution[5].Thickness - md.results.TransientSolution[5].Thickness,
                md2.results.TransientSolution[5].IceVolume - md.results.TransientSolution[5].IceVolume,
                md2.results.TransientSolution[5].Temperature - md.results.TransientSolution[5].Temperature,
                md2.results.TransientSolution[5].Pressure - md.results.TransientSolution[5].Pressure,
                md2.results.TransientSolution[6].Vx - md.results.TransientSolution[6].Vx,
                md2.results.TransientSolution[6].Vy - md.results.TransientSolution[6].Vy,
                md2.results.TransientSolution[6].Vel - md.results.TransientSolution[6].Vel,
                md2.results.TransientSolution[6].TotalSmb - md.results.TransientSolution[6].TotalSmb,
                md2.results.TransientSolution[6].Base - md.results.TransientSolution[6].Base,
                md2.results.TransientSolution[6].Surface - md.results.TransientSolution[6].Surface,
                md2.results.TransientSolution[6].Thickness - md.results.TransientSolution[6].Thickness,
                md2.results.TransientSolution[6].IceVolume - md.results.TransientSolution[6].IceVolume,
                md2.results.TransientSolution[6].Temperature - md.results.TransientSolution[6].Temperature,
                md2.results.TransientSolution[6].Pressure - md.results.TransientSolution[6].Pressure,
                md2.results.TransientSolution[7].Vx - md.results.TransientSolution[7].Vx,
                md2.results.TransientSolution[7].Vy - md.results.TransientSolution[7].Vy,
                md2.results.TransientSolution[7].Vel - md.results.TransientSolution[7].Vel,
                md2.results.TransientSolution[7].TotalSmb - md.results.TransientSolution[7].TotalSmb,
                md2.results.TransientSolution[7].Base - md.results.TransientSolution[7].Base,
                md2.results.TransientSolution[7].Surface - md.results.TransientSolution[7].Surface,
                md2.results.TransientSolution[7].Thickness - md.results.TransientSolution[7].Thickness,
                md2.results.TransientSolution[7].IceVolume - md.results.TransientSolution[7].IceVolume,
                md2.results.TransientSolution[7].Temperature - md.results.TransientSolution[7].Temperature,
                md2.results.TransientSolution[7].Pressure - md.results.TransientSolution[7].Pressure]
