#Test Name: 79NorthPISMhydrofriction
import pyissm
import numpy as np

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/79North.exp', 10000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/79NorthShelf.exp', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/79North.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3

# Hydrology
md.hydrology = pyissm.model.classes.hydrology.pism()
md.hydrology.drainage_rate = 0.001 * np.ones((md.mesh.numberofvertices))
md.hydrology.watercolumn_max = 2 * np.ones((md.mesh.numberofvertices))
md.initialization.pressure = np.zeros((md.mesh.numberofvertices))
md.initialization.watercolumn = np.zeros((md.mesh.numberofvertices))
md.basalforcings.groundedice_melting_rate = np.arange(0, md.mesh.numberofvertices) + 1

# Friction
md.friction = pyissm.model.classes.friction.pism()
md.friction.till_friction_angle = 30 * np.ones((md.mesh.numberofvertices))
md.friction.sediment_compressibility_coefficient = 0.12 * np.ones((md.mesh.numberofvertices))

md.transient.ishydrology = 1
md.transient.issmb = 0
md.transient.ismasstransport = 0
md.transient.isstressbalance = 1
md.transient.isthermal = 0
md.transient.isgroundingline = 0

# Execute model
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Vx1', 'Vx2', 'Vx3', 'Vy1', 'Vy2', 'Vy3', 'Vel1', 'Vel2', 'Vel3']
field_tolerances = [2e-11, 2e-11, 2e-11, 2e-11, 2e-11, 2e-11, 2e-11, 2e-11, 2e-11]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[1].Vel,
                md.results.TransientSolution[2].Vel]
