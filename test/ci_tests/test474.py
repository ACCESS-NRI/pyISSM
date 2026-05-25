#Test Name: PicoMeltRate_HO
import pyissm
import numpy as np

md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 90000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/SquareShelf.exp', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetShelf.py')
md.initialization.vx[:] = 1.
md.initialization.vy[:] = 1.
md.geometry.thickness[:] = 500. - md.mesh.x / 10000.
md.geometry.bed = -100. - md.mesh.x / 1000.
md.geometry.base = -md.geometry.thickness * md.materials.rho_ice / md.materials.rho_water
md.mask.ocean_levelset = md.geometry.thickness + md.materials.rho_water / md.materials.rho_ice * md.geometry.bed
pos = np.where(md.mask.ocean_levelset >= 0.)
md.geometry.base[pos] = md.geometry.bed[pos]
md.geometry.surface = md.geometry.base + md.geometry.thickness
md = md.extrude(3, 1.1)
md = pyissm.model.param.set_flow_equation(md, HO='all')

# Set PICO Parameters
md.basalforcings = pyissm.model.classes.basalforcings.pico()
md.basalforcings.basin_id = np.zeros((md.mesh.numberofelements,))
y_elem = np.mean(md.mesh.y[md.mesh.elements.astype(int) - 1], axis=1)
md.basalforcings.basin_id[y_elem >= 5e5] = 1
md.basalforcings.basin_id[y_elem < 5e5] = 2
md.basalforcings.num_basins = 2
md.basalforcings.farocean_temperature = np.array([[271.15, 272.15, 273.15], [274.15, 275.15, 276.15], [0.5, 1., 1.5]])
md.basalforcings.farocean_salinity = np.array([[31., 32., 33.], [34., 35., 36.], [0.5, 1., 1.5]])
md.basalforcings.maxboxcount = 5
md.basalforcings.isplume = 0

# Boundary conditions
md.mask.ice_levelset = -np.ones((md.mesh.numberofvertices,))
md.mask.ice_levelset[np.where(md.mesh.x == max(md.mesh.x))] = 0.

# Model conditions
md.transient.isthermal = 0
md.transient.isstressbalance = 1
md.transient.isgroundingline = 1
md.transient.ismasstransport = 1
md.transient.issmb = 1
md.transient.requested_outputs = ['default', 'BasalforcingsFloatingiceMeltingRate']
md.groundingline.migration = 'SubelementMigration'
md.groundingline.friction_interpolation = 'SubelementFriction1'
md.groundingline.melt_interpolation = 'SubelementMelt1'
md.timestepping.final_time = 1.5
md.timestepping.time_step = 0.5
md.cluster.np = 3

md = pyissm.model.execute.solve(md, 'Transient')

field_names = ['Bed1', 'Surface1', 'Thickness1', 'Floatingice1', 'Vx1', 'Vy1', 'Pressure1', 'FloatingiceMeltingrate1',
               'Bed2', 'Surface2', 'Thickness2', 'Floatingice2', 'Vx2', 'Vy2', 'Pressure2', 'FloatingiceMeltingrate2',
               'Bed3', 'Surface3', 'Thickness3', 'Floatingice3', 'Vx3', 'Vy3', 'Pressure3', 'FloatingiceMeltingrate3']
field_tolerances = [7e-9, 8e-9, 8e-9, 7e-9, 6e-8, 7e-8, 6e-9, 8e-10,
                    7e-9, 8e-9, 8e-9, 7e-9, 6e-8, 7e-8, 6e-9, 8e-10,
                    7e-9, 8e-9, 8e-9, 7e-9, 6e-8, 7e-8, 6e-9, 8e-10]
field_values = [md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].MaskOceanLevelset,
                md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].BasalforcingsFloatingiceMeltingRate,
                md.results.TransientSolution[1].Base,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].MaskOceanLevelset,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Pressure,
                md.results.TransientSolution[1].BasalforcingsFloatingiceMeltingRate,
                md.results.TransientSolution[2].Base,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].MaskOceanLevelset,
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Pressure,
                md.results.TransientSolution[2].BasalforcingsFloatingiceMeltingRate]
