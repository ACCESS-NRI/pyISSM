import inspect
import os
import numpy as np
import pyissm

#Geometry
hmin = 300.
hmax = 1000.
ymin = min(md.mesh.y)
ymax = max(md.mesh.y)
xmin = min(md.mesh.x)
xmax = max(md.mesh.x)
md.geometry.thickness = hmax + (hmin - hmax) * (md.mesh.y - ymin) / (ymax - ymin) + 0.1 * (hmin - hmax) * (md.mesh.x - xmin) / (xmax - xmin)
md.geometry.base = -md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness
md.geometry.surface = md.geometry.base + md.geometry.thickness

#Initial velocity and pressure
x = np.array(pyissm.tools.archive.arch_read('../assets/Data/SquareShelf.arch', 'x'))
y = np.array(pyissm.tools.archive.arch_read('../assets/Data/SquareShelf.arch', 'y'))
vx = np.array(pyissm.tools.archive.arch_read('../assets/Data/SquareShelf.arch', 'vx'))
vy = np.array(pyissm.tools.archive.arch_read('../assets/Data/SquareShelf.arch', 'vy'))
index = pyissm.tools.archive.arch_read('../assets/Data/SquareShelf.arch', 'index').astype(int)
# Initial velocity and pressure
md.initialization.vx = pyissm.tools.wrappers.InterpFromMeshToMesh2d(index, x, y, vx, md.mesh.x, md.mesh.y)
md.initialization.vy = pyissm.tools.wrappers.InterpFromMeshToMesh2d(index, x, y, vy, md.mesh.x, md.mesh.y)
md.initialization.vz = np.zeros((md.mesh.numberofvertices))
md.initialization.pressure = np.zeros((md.mesh.numberofvertices))

# Materials
md.initialization.temperature = (273. - 20.) * np.ones((md.mesh.numberofvertices))
md.materials.rheology_B = pyissm.tools.materials.paterson(md.initialization.temperature)
md.materials.rheology_n = 3. * np.ones((md.mesh.numberofelements))

#Friction
md.friction.coefficient = 20. * np.ones((md.mesh.numberofvertices))
md.friction.coefficient[np.nonzero(md.mask.ocean_levelset < 0.)[0]] = 0.
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))

# Numerical parameters
md.masstransport.stabilization = 1.
md.thermal.stabilization = 1.
md.settings.waitonlock = 30
md.stressbalance.restol = 0.10
md.steadystate.reltol = 0.02
md.stressbalance.reltol = 0.02
md.stressbalance.abstol = np.nan
md.timestepping.time_step = 1.
md.timestepping.final_time = 3.
md.groundingline.migration = 'None'

# Boundary conditions:
md = pyissm.model.bc.set_ice_shelf_bc(md, '../assets/Exp/SquareFront2.exp')

# Change name so that no test have the same name
if len(inspect.stack()) > 2:
    md.miscellaneous.name = os.path.basename(inspect.stack()[2][1]).split('.')[0]
