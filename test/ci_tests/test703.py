#Test Name: FlowbandFSsheetshelfTrans
import pyissm
import numpy as np
import copy

# mesh parameters
x = np.arange(-5, 5.5, .5).T
[b, h, sea] = pyissm.tools.geometry.nowicki_profile(x)
x = x * 10**3
h = h * 10**3
b = (b - sea) * 10**3

# mesh domain
md = pyissm.model.mesh.bamg_flowband(pyissm.model.Model(), x, b + h, b, hmax = 150.)

# geometry
md.geometry.surface = np.interp(md.mesh.x, x, b + h)
md.geometry.base = np.interp(md.mesh.x, x, b)
md.geometry.thickness = md.geometry.surface - md.geometry.base

# mask
md.mask.ice_levelset = -np.ones((md.mesh.numberofvertices))
md.mask.ice_levelset[np.where(md.mesh.vertex_flags(2))] = 0
md.mask.ocean_levelset = -0.5 * np.ones((md.mesh.numberofvertices))
md.mask.ocean_levelset[np.where(md.mesh.x < 0)] = 0.5

# materials
md.initialization.temperature = (273. - 20.) * np.ones((md.mesh.numberofvertices))
md.materials.rheology_B = pyissm.tools.materials.paterson(md.initialization.temperature)
md.materials.rheology_n = 3 * np.ones((md.mesh.numberofelements))

# friction
md.friction.coefficient = np.zeros((md.mesh.numberofvertices))
md.friction.coefficient[np.where(md.mesh.vertex_flags(1))] = 20
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))

# boundary conditions
md.stressbalance.spcvx = np.nan * np.ones((md.mesh.numberofvertices))
md.stressbalance.spcvy = np.nan * np.ones((md.mesh.numberofvertices))
md.stressbalance.spcvz = np.nan * np.ones((md.mesh.numberofvertices))
md.stressbalance.referential = np.nan * np.ones((md.mesh.numberofvertices, 6))
md.stressbalance.loadingforce = np.zeros((md.mesh.numberofvertices, 3))
md.stressbalance.spcvx[np.where(md.mesh.vertex_flags(4))] = 800.
md.stressbalance.spcvy[np.where(md.mesh.vertex_flags(4))] = 0.

# Misc
md = pyissm.model.param.set_flow_equation(md, FS = 'all')
md.flowequation.fe_FS = 'TaylorHood'
md.stressbalance.abstol = np.nan
md.miscellaneous.name = 'test703'
md.cluster.np = 3

# Transient settings
md.timestepping.time_step = 0.000001
md.timestepping.final_time = 0.000005
md.stressbalance.shelf_dampening = 1.
md.smb.mass_balance = np.zeros((md.mesh.numberofvertices))
md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices))
md.basalforcings.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices))
md.basalforcings.geothermalflux = np.zeros((md.mesh.numberofvertices))
posb = np.where(np.logical_and(md.mesh.x > 0., md.mesh.vertexonbase))[0]
md.basalforcings.groundedice_melting_rate[posb] = 18.
md.basalforcings.floatingice_melting_rate[posb] = 18.
md.initialization.vx = np.zeros((md.mesh.numberofvertices))
md.initialization.vy = np.zeros((md.mesh.numberofvertices))
md.initialization.pressure = np.zeros((md.mesh.numberofvertices))
md.masstransport.spcthickness = np.nan * np.ones((md.mesh.numberofvertices))
md.thermal.spctemperature = np.nan * np.ones((md.mesh.numberofvertices))
md.transient.isthermal = 0
md.masstransport.isfreesurface = 1
md.groundingline.migration = 'None'

# Execute model(s)
md.stressbalance.shelf_dampening = 1
md1 = copy.deepcopy(md)
md1 = pyissm.model.execute.solve(md1, 'Transient')

md.stressbalance.shelf_dampening = 0
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'Pressure1', 'Bed1', 'Surface1', 'Thickness1',
               'Vx2', 'Vy2', 'Vel2', 'Pressure2', 'Bed2', 'Surface2', 'Thickness2',
               'Vx3', 'Vy3', 'Vel3', 'Pressure3', 'Bed3', 'Surface3', 'Thickness3',
               'Vx1_damp', 'Vy1_damp', 'Vel1_damp', 'Pressure1_damp', 'Bed1_damp', 'Surface1_damp', 'Thickness1_damp',
               'Vx2_damp', 'Vy2_damp', 'Vel2_damp', 'Pressure2_damp', 'Bed2_damp', 'Surface2_damp', 'Thickness2_damp',
               'Vx3_damp', 'Vy3_damp', 'Vel3_damp', 'Pressure3_damp', 'Bed3_damp', 'Surface3_damp', 'Thickness3_damp']
field_tolerances = [2e-08, 2e-08, 2e-08, 1e-08, 1e-10, 1e-10, 1e-10,
                    2e-08, 2e-08, 2e-08, 1e-08, 1e-10, 1e-10, 1e-10,
                    2e-08, 2e-08, 2e-08, 1e-08, 1e-10, 1e-10, 1e-10,
                    5e-08, 5e-08, 5e-08, 1e-08, 1e-10, 1e-10, 1e-10,
                    5e-08, 5e-08, 5e-08, 1e-08, 1e-10, 1e-10, 1e-10,
                    5e-08, 5e-08, 5e-08, 1e-08, 1e-10, 1e-10, 1e-10]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Vel,
                md.results.TransientSolution[1].Pressure,
                md.results.TransientSolution[1].Base,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Vel,
                md.results.TransientSolution[2].Pressure,
                md.results.TransientSolution[2].Base,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].Thickness,
                md1.results.TransientSolution[0].Vx,
                md1.results.TransientSolution[0].Vy,
                md1.results.TransientSolution[0].Vel,
                md1.results.TransientSolution[0].Pressure,
                md1.results.TransientSolution[0].Base,
                md1.results.TransientSolution[0].Surface,
                md1.results.TransientSolution[0].Thickness,
                md1.results.TransientSolution[1].Vx,
                md1.results.TransientSolution[1].Vy,
                md1.results.TransientSolution[1].Vel,
                md1.results.TransientSolution[1].Pressure,
                md1.results.TransientSolution[1].Base,
                md1.results.TransientSolution[1].Surface,
                md1.results.TransientSolution[1].Thickness,
                md1.results.TransientSolution[2].Vx,
                md1.results.TransientSolution[2].Vy,
                md1.results.TransientSolution[2].Vel,
                md1.results.TransientSolution[2].Pressure,
                md1.results.TransientSolution[2].Base,
                md1.results.TransientSolution[2].Surface,
                md1.results.TransientSolution[2].Thickness]
