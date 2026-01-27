import os.path
import inspect
import pyissm
import numpy as np

# Geometry and observation
x = np.array(pyissm.tools.archive.arch_read('../assets/Data/79North.arch', 'x'))
y = np.array(pyissm.tools.archive.arch_read('../assets/Data/79North.arch', 'y'))
vx = np.array(pyissm.tools.archive.arch_read('../assets/Data/79North.arch', 'vx'))
vy = np.array(pyissm.tools.archive.arch_read('../assets/Data/79North.arch', 'vy'))
index = np.array(pyissm.tools.archive.arch_read('../assets/Data/79North.arch', 'index')).astype(int)
surface = np.array(pyissm.tools.archive.arch_read('../assets/Data/79North.arch', 'surface'))
thickness = np.array(pyissm.tools.archive.arch_read('../assets/Data/79North.arch', 'thickness'))

md.initialization.vx = pyissm.tools.wrappers.InterpFromMeshToMesh2d(index, x, y, vx, md.mesh.x, md.mesh.y)
md.initialization.vy = pyissm.tools.wrappers.InterpFromMeshToMesh2d(index, x, y, vy, md.mesh.x, md.mesh.y)
md.geometry.surface = pyissm.tools.wrappers.InterpFromMeshToMesh2d(index, x, y, surface, md.mesh.x, md.mesh.y)
md.geometry.thickness = pyissm.tools.wrappers.InterpFromMeshToMesh2d(index, x, y, thickness, md.mesh.x, md.mesh.y)
md.geometry.base = md.geometry.surface - md.geometry.thickness

# Materials
md.initialization.temperature = (273. - 20.) * np.ones((md.mesh.numberofvertices))
md.materials.rheology_B = pyissm.tools.materials.paterson(md.initialization.temperature)
md.materials.rheology_n = 3. * np.ones((md.mesh.numberofelements))
md.initialization.temperature = md.initialization.temperature

# Friction
md.friction.coefficient = 50. * np.ones((md.mesh.numberofvertices))
md.friction.coefficient[np.nonzero(md.mask.ocean_levelset < 0.)[0]] = 0.
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))

# Ice shelf melting and surface mass balance
md.basalforcings.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices))
md.basalforcings.floatingice_melting_rate[np.nonzero(md.mask.ocean_levelset < 0.)[0]] = 0.
md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices))
md.smb.mass_balance = 15 * np.ones((md.mesh.numberofvertices))

# Numerical parameters
md.masstransport.stabilization = 1
md.thermal.stabilization = 1
md.settings.waitonlock = 30
md.timestepping.time_step = 1.
md.timestepping.final_time = 3.
md.stressbalance.restol = 0.05
md.stressbalance.reltol = 0.005
md.steadystate.reltol = 0.005
md.stressbalance.abstol = np.nan
md.groundingline.migration = 'None'

# Boundary conditions:
md = pyissm.model.bc.set_marine_ice_sheet_bc(md)
pos = np.nonzero(md.mesh.vertexonboundary)
md.balancethickness.spcthickness[pos] = md.geometry.thickness[pos]
md.masstransport.spcthickness[pos] = md.geometry.thickness[pos]

# Change name so that no test have the same name
if len(inspect.stack()) > 2:
    md.miscellaneous.name = os.path.basename(inspect.stack()[2][1]).split('.')[0]
