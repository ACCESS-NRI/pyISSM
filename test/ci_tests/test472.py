#Test Name: ISMIP6MeltRateTest
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
md = pyissm.model.param.set_flow_equation(md, SSA='all')

# Set ISMIP6 melt rate parameters
md.basalforcings = pyissm.model.classes.basalforcings.ismip6()
md.basalforcings.basin_id = np.zeros((md.mesh.numberofelements,))
y_elem = np.mean(md.mesh.y[md.mesh.elements.astype(int) - 1], axis=1)
md.basalforcings.basin_id[y_elem >= 5e5] = 1
md.basalforcings.basin_id[y_elem < 5e5] = 2
md.basalforcings.num_basins = 2
md.basalforcings.delta_t = np.array([0.1, 0.2])
md.basalforcings.tf_depths = np.array([0., -1000., -2000.])
md.basalforcings.gamma_0 = 14477.
md.basalforcings.islocal = 0

# Build an artificial tf field (for times 0 and 1, 3 depth layers)
nv = md.mesh.numberofvertices
A = np.vstack([np.column_stack([1. * np.ones(nv), 1.5 * np.ones(nv)]), [0., 1.]])
B = np.vstack([np.column_stack([2. * np.ones(nv), 2.5 * np.ones(nv)]), [0., 1.]])
C = np.vstack([np.column_stack([3. * np.ones(nv), 3.5 * np.ones(nv)]), [0., 1.]])
md.basalforcings.tf = [A, B, C]

# Add a melt anomaly
md.basalforcings.melt_anomaly = np.vstack([np.column_stack([np.ones(nv), 2. * np.ones(nv)]), [1., 2.]])

# Boundary conditions
md.mask.ice_levelset = -np.ones((md.mesh.numberofvertices,))
md.mask.ice_levelset[np.where(md.mesh.x == max(md.mesh.x))] = 0.

# Model conditions
md.transient.isthermal = 0
md.transient.isstressbalance = 1
md.transient.isgroundingline = 1
md.transient.ismasstransport = 1
md.transient.issmb = 1
md.transient.requested_outputs = ['default', 'BasalforcingsFloatingiceMeltingRate', 'BasalforcingsIsmip6TfShelf']
md.groundingline.migration = 'SubelementMigration'
md.groundingline.friction_interpolation = 'SubelementFriction1'
md.groundingline.melt_interpolation = 'SubelementMelt1'
md.timestepping.final_time = 1.5
md.timestepping.time_step = 0.5
md.cluster.np = 1

md.basalforcings.initialize(md)
md = pyissm.model.execute.solve(md, 'Transient')

field_names = ['Bed1', 'Surface1', 'Thickness1', 'Floatingice1', 'Vx1', 'Vy1', 'Pressure1', 'FloatingiceMeltingrate1', 'ThermalForcing1',
               'Bed2', 'Surface2', 'Thickness2', 'Floatingice2', 'Vx2', 'Vy2', 'Pressure2', 'FloatingiceMeltingrate2', 'ThermalForcing2',
               'Bed3', 'Surface3', 'Thickness3', 'Floatingice3', 'Vx3', 'Vy3', 'Pressure3', 'FloatingiceMeltingrate3', 'ThermalForcing3']
field_tolerances = [7e-9, 8e-9, 8e-9, 7e-9, 6e-8, 7e-8, 6e-9, 8e-10, 7e-8,
                    7e-9, 8e-9, 8e-9, 7e-9, 6e-8, 7e-8, 6e-9, 8e-10, 7e-8,
                    7e-9, 8e-9, 8e-9, 7e-9, 6e-8, 7e-8, 6e-9, 8e-10, 7e-8]
field_values = [md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].MaskOceanLevelset,
                md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].BasalforcingsFloatingiceMeltingRate,
                md.results.TransientSolution[0].BasalforcingsIsmip6TfShelf,
                md.results.TransientSolution[1].Base,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].MaskOceanLevelset,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Pressure,
                md.results.TransientSolution[1].BasalforcingsFloatingiceMeltingRate,
                md.results.TransientSolution[1].BasalforcingsIsmip6TfShelf,
                md.results.TransientSolution[2].Base,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].MaskOceanLevelset,
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Pressure,
                md.results.TransientSolution[2].BasalforcingsFloatingiceMeltingRate,
                md.results.TransientSolution[2].BasalforcingsIsmip6TfShelf]
