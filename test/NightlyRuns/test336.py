#Test Name: SquareSheetConstrainedSmbComponents2d
import pyissm
import numpy as np

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000.)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelf.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3
md.timestepping.time_step = 1.
md.settings.output_frequency = 1
md.timestepping.final_time = 4.

# Set up transient
smb = np.ones((md.mesh.numberofvertices)) * 3.6
smb = np.vstack((smb, smb * -1.)).T

md.smb = pyissm.model.classes.smb.components()
md.smb.accumulation = np.vstack((smb * 2, [1.5, 3.]))
md.smb.runoff = np.vstack((smb / 2, [1.5, 3.]))
md.smb.evaporation = np.vstack((smb / 2, [1.5, 3.]))
md.transient.isthermal = False

# Execute model
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'Pressure1', 'Bed1', 'Surface1', 'Thickness1', 'SmbMassBalance1',
               'Vx2', 'Vy2', 'Vel2', 'Pressure2', 'Bed2', 'Surface2', 'Thickness2', 'SmbMassBalance2',
               'Vx3', 'Vy3', 'Vel3', 'Pressure3', 'Bed3', 'Surface3', 'Thickness3', 'SmbMassBalance3',
               'Vx4', 'Vy4', 'Vel4', 'Pressure4', 'Bed4', 'Surface4', 'Thickness4', 'SmbMassBalance4']
field_tolerances = [1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10,
                    1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10,
                    1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10,
                    1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].SmbMassBalance,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Vel,
                md.results.TransientSolution[1].Pressure,
                md.results.TransientSolution[1].Base,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].SmbMassBalance,
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Vel,
                md.results.TransientSolution[2].Pressure,
                md.results.TransientSolution[2].Base,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].SmbMassBalance,
                md.results.TransientSolution[3].Vx,
                md.results.TransientSolution[3].Vy,
                md.results.TransientSolution[3].Vel,
                md.results.TransientSolution[3].Pressure,
                md.results.TransientSolution[3].Base,
                md.results.TransientSolution[3].Surface,
                md.results.TransientSolution[3].Thickness,
                md.results.TransientSolution[3].SmbMassBalance]
