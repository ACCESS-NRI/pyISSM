#Test Name: FlowbandFSsheetshelf
import pyissm
import numpy as np

# mesh parameters
x = np.arange(-5, 5.5, .5).T
[b, h, sea] = pyissm.tools.geometry.nowicki_profile(x)
x = x * 10**3
h = h * 10**3
b = (b - sea) * 10**3

# mesh domain
md = pyissm.model.mesh.bamg_flowband(pyissm.model.Model(), x, b + h, b, hmax = 150.)

# parameterize
md.geometry.surface = np.interp(md.mesh.x, x, b + h)
md.geometry.base = np.interp(md.mesh.x, x, b)
md.geometry.thickness = md.geometry.surface - md.geometry.base

# mask
md.mask.ice_levelset = -np.ones((md.mesh.numberofvertices, ))
md.mask.ice_levelset[np.where(md.mesh.vertex_flags(2))] = 0
md.mask.ocean_levelset = -0.5 * np.ones((md.mesh.numberofvertices))
md.mask.ocean_levelset[np.where(md.mesh.x < 0)] = 0.5

# materials
md.initialization.temperature = (273. - 20.) * np.ones((md.mesh.numberofvertices, ))
md.materials.rheology_B = pyissm.tools.materials.paterson(md.initialization.temperature)
md.materials.rheology_n = 3 * np.ones((md.mesh.numberofelements, ))

# damage
md.damage.D = np.zeros((md.mesh.numberofvertices, ))
md.damage.spcdamage = float('NaN') * np.ones((md.mesh.numberofvertices, ))

# friciton
md.friction.coefficient = np.zeros((md.mesh.numberofvertices, ))
md.friction.coefficient[np.where(md.mesh.vertex_flags(1))] = 20
md.friction.p = np.ones((md.mesh.numberofelements, ))
md.friction.q = np.ones((md.mesh.numberofelements, ))

# boundary conditions
md.stressbalance.spcvx = np.nan * np.ones((md.mesh.numberofvertices, ))
md.stressbalance.spcvy = np.nan * np.ones((md.mesh.numberofvertices, ))
md.stressbalance.spcvz = np.nan * np.ones((md.mesh.numberofvertices, ))
md.stressbalance.referential = np.nan * np.ones((md.mesh.numberofvertices, 6))
md.stressbalance.loadingforce = np.zeros((md.mesh.numberofvertices, 3))
md.stressbalance.spcvx[np.where(md.mesh.vertex_flags(4))] = 800.
md.stressbalance.spcvy[np.where(md.mesh.vertex_flags(4))] = 0.
md.basalforcings.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices, ))

# misc
md = pyissm.model.param.set_flow_equation(md, FS = 'all')
md.stressbalance.abstol = np.nan
md.stressbalance.FSreconditioning = 1
md.stressbalance.maxiter = 20
md.flowequation.augmented_lagrangian_r = 10000.
md.flowequation.augmented_lagrangian_rhop = 10000.
md.initialization.pressure = md.constants.g * md.materials.rho_ice * (md.geometry.surface - md.mesh.y)
md.miscellaneous.name = 'test702'
md.groundingline.migration = 'None'
md.cluster.np = 2

# Execute model
field_names = []
field_tolerances = []
field_values = []
for i in ['MINI', 'MINIcondensed', 'TaylorHood', 'XTaylorHood', 'LATaylorHood']:
    print(' ')
    print(' == == == Testing ' + i + ' Full - Stokes Finite element == == = ')
    md.flowequation.fe_FS = i
    md = pyissm.model.execute.solve(md, 'Stressbalance')
    field_names.extend(['Vx' + i, 'Vy' + i, 'Vel' + i, 'Pressure' + i])
    field_tolerances.extend([8e-5, 8e-5, 8e-5, 1e-08])
    field_values.extend([md.results.StressbalanceSolution.Vx,
                         md.results.StressbalanceSolution.Vy,
                         md.results.StressbalanceSolution.Vel,
                         md.results.StressbalanceSolution.Pressure])
