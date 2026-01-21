<<<<<<< HEAD

=======
>>>>>>> 184d2d6 (Add test101; Update massfluxatgate & outputdef...)
import os.path
import pyissm
import numpy as np
import inspect

<<<<<<< HEAD

#Start defining model parameters here
#Geometry
=======
# Geometry
>>>>>>> 184d2d6 (Add test101; Update massfluxatgate & outputdef...)
hmin = 300.
hmax = 1000.
ymin = np.min(md.mesh.y)
ymax = np.max(md.mesh.y)
xmin = min(md.mesh.x)
xmax = max(md.mesh.x)
md.geometry.thickness = hmax + (hmin - hmax) * (md.mesh.y - ymin) / (ymax - ymin) + 0.1 * (hmin - hmax) * (md.mesh.x - xmin) / (xmax - xmin)
md.geometry.base = -md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness
md.geometry.surface = md.geometry.base + md.geometry.thickness
md.geometry.bed = md.geometry.base - 10

<<<<<<< HEAD
#Initial velocity
=======
# Initial velocity
>>>>>>> 184d2d6 (Add test101; Update massfluxatgate & outputdef...)
x = np.array(pyissm.tools.archive.arch_read('../assets/Data/SquareShelfConstrained.arch', 'x'))
y = np.array(pyissm.tools.archive.arch_read('../assets/Data/SquareShelfConstrained.arch', 'y'))
vx = np.array(pyissm.tools.archive.arch_read('../assets/Data/SquareShelfConstrained.arch', 'vx'))
vy = np.array(pyissm.tools.archive.arch_read('../assets/Data/SquareShelfConstrained.arch', 'vy'))
index = np.array(pyissm.tools.archive.arch_read('../assets/Data/SquareShelfConstrained.arch', 'index').astype(int))
md.initialization.vx = pyissm.tools.wrappers.InterpFromMeshToMesh2d(index, x, y, vx, md.mesh.x, md.mesh.y)
md.initialization.vy = pyissm.tools.wrappers.InterpFromMeshToMesh2d(index, x, y, vy, md.mesh.x, md.mesh.y)
md.initialization.vz = np.zeros((md.mesh.numberofvertices))
md.initialization.pressure = np.zeros((md.mesh.numberofvertices))

<<<<<<< HEAD
#Materials
=======
# Materials
>>>>>>> 184d2d6 (Add test101; Update massfluxatgate & outputdef...)
md.initialization.temperature = (273. - 20.) * np.ones((md.mesh.numberofvertices))
md.materials.rheology_B = pyissm.tools.materials.paterson(md.initialization.temperature)
md.materials.rheology_n = 3. * np.ones((md.mesh.numberofelements))

<<<<<<< HEAD
#Surface mass balance and basal melting
=======
# Surface mass balance and basal melting
>>>>>>> 184d2d6 (Add test101; Update massfluxatgate & outputdef...)
md.smb.mass_balance = 10. * np.ones((md.mesh.numberofvertices))
md.basalforcings.groundedice_melting_rate = 5. * np.ones((md.mesh.numberofvertices))
md.basalforcings.floatingice_melting_rate = 5. * np.ones((md.mesh.numberofvertices))

<<<<<<< HEAD
#Friction
=======
# Friction
>>>>>>> 184d2d6 (Add test101; Update massfluxatgate & outputdef...)
md.friction.coefficient = 20. * np.ones((md.mesh.numberofvertices))
md.friction.coefficient[np.nonzero(md.mask.ocean_levelset < 0.)[0]] = 0.
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))

<<<<<<< HEAD
#Numerical parameters
=======
# Numerical parameters
>>>>>>> 184d2d6 (Add test101; Update massfluxatgate & outputdef...)
md.masstransport.stabilization = 1
md.thermal.stabilization = 1
md.settings.waitonlock = 30
md.stressbalance.restol = 0.05
md.stressbalance.reltol = 0.05
md.steadystate.reltol = 0.05
md.stressbalance.abstol = float('nan')
md.timestepping.time_step = 1.
md.timestepping.final_time = 3.

<<<<<<< HEAD
#Deal with boundary conditions:
md = pyissm.model.bc.set_ice_shelf_bc(md)

#Change name so that no tests have the same name
if len(inspect.stack()) > 2:
    md.miscellaneous.name = os.path.basename(inspect.stack()[2][1]).split('.')[0]
=======
# Deal with boundary conditions:
md = pyissm.model.bc.set_ice_shelf_bc(md)

# Change name so that no tests have the same name
if len(inspect.stack()) > 2:
    md.miscellaneous.name = os.path.basename(inspect.stack()[2][1]).split('.')[0]
>>>>>>> 184d2d6 (Add test101; Update massfluxatgate & outputdef...)
