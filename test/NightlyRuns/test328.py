#Test Name: SquareSheetConstrainedSmbGradients2d
import pyissm
import copy

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000.)
md = pyissm.model.param.set_mask(md, None, None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetConstrained.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3

md.smb = pyissm.model.classes.smb.gradients()
md.smb.b_pos = (-100. + 0.00005 * md.mesh.x - 0.0001 * md.mesh.y) / 1000. * md.materials.rho_freshwater / md.materials.rho_ice
md.smb.b_neg = (250. + 0.000051 * md.mesh.x - 0.00011 * md.mesh.y) / 1000. * md.materials.rho_freshwater / md.materials.rho_ice
md.transient.requested_outputs = ['default', 'TotalSmb']
md.smb.href = copy.deepcopy(md.geometry.surface)
md.smb.smbref= (1000. - 0.001 * md.mesh.x - 0.005 * md.mesh.y) / 1000. * md.materials.rho_freshwater / md.materials.rho_ice

# Execute model
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1',
               'Bed1', 'Surface1', 'Thickness1',
               'SMB1', 'TotalSmb1',
               'Vx2', 'Vy2', 'Vel2', 'Bed2',
               'Surface2', 'Thickness2',
               'SMB2', 'TotalSmb2', 'Vx3', 'Vy3',
               'Vel3', 'Bed3', 'Surface3',
               'Thickness3', 'SMB3', 'TotalSmb3']
field_tolerances = [1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13,
                    1e-13, 2e-13, 1e-13]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].SmbMassBalance,
                md.results.TransientSolution[0].TotalSmb,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Vel,
                md.results.TransientSolution[1].Base,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].TotalSmb,
                md.results.TransientSolution[1].SmbMassBalance,
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Vel,
                md.results.TransientSolution[2].Base,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].SmbMassBalance,
                md.results.TransientSolution[2].TotalSmb]
