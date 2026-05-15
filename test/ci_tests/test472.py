#Test Name: ISMIP6MeltRateTest 
import pyissm
import numpy as np

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 90000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/SquareShelf.exp', '')
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetShelf.par')

md.initialization.vx[:] = 1.
md.initialization.vy[:] = 1.
md.geometry.thickness[:] = 500. - md.mesh.x / 10000.
md.geometry.bed = -100. - md.mesh.x / 1000.
md.geometry.base = (-md.geometry.thickness * md.materials.rho_ice / md.materials.rho_water)
md.mask.ocean_levelset = (md.geometry.thickness + md.materials.rho_water / md.materials.rho_ice * md.geometry.bed)
pos = np.where(md.mask.ocean_levelset >= 0.)[0]
md.geometry.base[pos] = md.geometry.bed[pos]
md.geometry.surface = (md.geometry.base + md.geometry.thickness)

md = pyissm.model.param.set_flow_equation(md, SSA = 'all')

# ISMIP6 basal forcings
md.basalforcings = pyissm.model.classes.basalforcings.ismip6(md.basalforcings)
md.basalforcings.basin_id = np.zeros((md.mesh.numberofelements,))
yE = np.mean(md.mesh.y[md.mesh.elements], axis = 1)
pos1 = np.where(yE >= 5e5)[0]
md.basalforcings.basin_id[pos1] = 1
pos2 = np.where(yE < 5e5)[0]
md.basalforcings.basin_id[pos2] = 2
md.basalforcings.num_basins = 2
md.basalforcings.delta_t = np.array([0.1, 0.2])
md.basalforcings.tf_depths = np.array([0., -1000., -2000.])
md.basalforcings.gamma_0 = 14477.
md.basalforcings.islocal = 0


# Thermal forcing fields
# Time series:
# rows = vertices + optional time row
# cols = timesteps

temp1a = np.ones((md.mesh.numberofvertices,))
temp1b = 1.5 * np.ones((md.mesh.numberofvertices,))
A = np.vstack([np.column_stack([temp1a, temp1b]), np.array([[0., 1.]])])
temp2a = 2.0 * np.ones((md.mesh.numberofvertices,))
temp2b = 2.5 * np.ones((md.mesh.numberofvertices,))
B = np.vstack([np.column_stack([temp2a, temp2b]), np.array([[0., 1.]])])
temp3a = 3.0 * np.ones((md.mesh.numberofvertices,))
temp3b = 3.5 * np.ones((md.mesh.numberofvertices,))
C = np.vstack([np.column_stack([temp3a, temp3b]), np.array([[0., 1.]])])
md.basalforcings.tf = [A, B, C]

# Melt anomaly
md.basalforcings.melt_anomaly = np.ones((md.mesh.numberofvertices + 1, 2))
md.basalforcings.melt_anomaly[:, 1] = 2.
md.basalforcings.melt_anomaly[-1, :] = [1., 2.]

# Boundary conditions
md.mask.ice_levelset = -np.ones((md.mesh.numberofvertices, ))
pos = np.where(md.mesh.x == np.max(md.mesh.x))[0]
md.mask.ice_levelset[pos] = 0.

# Model conditions
md.transient.isthermal = 0
md.transient.isstressbalance = 1
md.transient.isgroundingline = 1
md.transient.ismasstransport = 1
md.transient.issmb = 1

md.transient.requested_outputs = [
    'default',
    'BasalforcingsFloatingiceMeltingRate',
    'BasalforcingsIsmip6TfShelf'
]

md.groundingline.migration = 'SubelementMigration'
md.groundingline.friction_interpolation = 'SubelementFriction1'
md.groundingline.melt_interpolation = 'SubelementMelt1'
md.timestepping.final_time = 1.5
md.timestepping.time_step = 0.5


# Solve
md.cluster.np = 1
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances
field_names = [
    'Bed1',
    'Surface1',
    'Thickness1',
    'Floatingice1',
    'Vx1',
    'Vy1',
    'Pressure1',
    'FloatingiceMeltingrate1',
    'ThermalForcing1',
    'Bed2',
    'Surface2',
    'Thickness2',
    'Floatingice2',
    'Vx2',
    'Vy2',
    'Pressure2',
    'FloatingiceMeltingrate2',
    'ThermalForcing2',
    'Bed3',
    'Surface3',
    'Thickness3',
    'Floatingice3',
    'Vx3',
    'Vy3',
    'Pressure3',
    'FloatingiceMeltingrate3',
    'ThermalForcing3',
]

field_tolerances = [
    7e-09, 8e-09, 8e-09, 7e-09,
    6e-08, 7e-08, 6e-09, 8e-10, 7e-08,
    7e-09, 8e-09, 8e-09, 7e-09,
    6e-08, 7e-08, 6e-09, 8e-10, 7e-08,
    7e-09, 8e-09, 8e-09, 7e-09,
    6e-08, 7e-08, 6e-09, 8e-10, 7e-08,
]

sol = md.results.TransientSolution
field_values = [
    sol[0].Base,
    sol[0].Surface,
    sol[0].Thickness,
    sol[0].MaskOceanLevelset,
    sol[0].Vx,
    sol[0].Vy,
    sol[0].Pressure,
    sol[0].BasalforcingsFloatingiceMeltingRate,
    sol[0].BasalforcingsIsmip6TfShelf,
    sol[1].Base,
    sol[1].Surface,
    sol[1].Thickness,
    sol[1].MaskOceanLevelset,
    sol[1].Vx,
    sol[1].Vy,
    sol[1].Pressure,
    sol[1].BasalforcingsFloatingiceMeltingRate,
    sol[1].BasalforcingsIsmip6TfShelf,
    sol[2].Base,
    sol[2].Surface,
    sol[2].Thickness,
    sol[2].MaskOceanLevelset,
    sol[2].Vx,
    sol[2].Vy,
    sol[2].Pressure,
    sol[2].BasalforcingsFloatingiceMeltingRate,
    sol[2].BasalforcingsIsmip6TfShelf,
]