#Test Name: SquareSheetConstrainedSmbGradientsEla2d
import pyissm
import numpy as np

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000.)
md = pyissm.model.param.set_mask(md, None, None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetConstrained.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')

md.smb = pyissm.model.classes.smb.gradientsela()
md.smb.ela = 1500. * np.ones((md.mesh.numberofvertices + 1, ))
md.smb.b_pos = 0.002 * np.ones((md.mesh.numberofvertices + 1, ))
md.smb.b_neg = 0.005 * np.ones((md.mesh.numberofvertices + 1, ))
md.smb.b_max = 4. * (md.materials.rho_freshwater / md.materials.rho_ice) * np.ones((md.mesh.numberofvertices + 1, ))
md.smb.b_min = -4. * (md.materials.rho_freshwater / md.materials.rho_ice) * np.ones((md.mesh.numberofvertices + 1, ))

# Change geometry
md.geometry.thickness = md.geometry.surface * 30.
md.geometry.surface = md.geometry.base + md.geometry.thickness

# Transient options
md.transient.requested_outputs = ['default', 'TotalSmb']

# Execute model
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'Bed1', 'Surface1', 'Thickness1', 'SMB1', 'TotalSmb1',
               'Vx2', 'Vy2', 'Vel2', 'Bed2', 'Surface2', 'Thickness2', 'SMB2', 'TotalSmb2',
               'Vx3', 'Vy3', 'Vel3', 'Bed3', 'Surface3', 'Thickness3', 'SMB3', 'TotalSmb3']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13,
                    1e-12, 1e-12, 1e-12, 1e-13, 1e-13, 1e-13, 1.5e-13, 1e-13]
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
