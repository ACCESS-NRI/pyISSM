#Test Name: SquareNoDynHydrologyDCOneLayer
import pyissm
import numpy as np

md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 100000.)
md = pyissm.model.param.set_mask(md, None, None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareNoDyn.py')
md.cluster.np = 1

md.transient.ishydrology = True
md.hydrology = pyissm.model.classes.hydrology.dc()
md.hydrology = md.hydrology.initialize(md)

md.hydrology.isefficientlayer = 0
md.hydrology.mask_thawed_node = np.ones((md.mesh.numberofvertices))
md.hydrology.sedimentlimit_flag = 1
md.hydrology.sedimentlimit = 8000.0

md.initialization.sediment_head = np.zeros((md.mesh.numberofvertices))
md.hydrology.spcsediment_head = np.nan * np.ones((md.mesh.numberofvertices))
pos = np.nonzero(md.mesh.y == 0.)[0]
md.hydrology.spcsediment_head[pos] = 0.0
md.basalforcings.groundedice_melting_rate = 2.0 * np.ones((md.mesh.numberofvertices))
md.basalforcings.floatingice_melting_rate = 0.0 * np.ones((md.mesh.numberofvertices))
md.hydrology.sediment_transmitivity = 3.0 * np.ones((md.mesh.numberofvertices))
md.timestepping.time_step = 0
md.timestepping.final_time = 1.0

# Execute model
md = pyissm.model.execute.solve(md, 'Hydrology')

#Fields and tolerances to track changes
#you can also compare with an analitic solution, but it is exact
#only if no limits are applied
#analitic=(md.mesh.y**2 - 2 * md.mesh.y * 1.0e6) * (-2.0 / (2 * md.constants.yts * md.hydrology.sediment_transmitivity))
field_names = ['SedimentWaterHead', 'SedimentHeadResidual']
field_tolerances = [1e-13, 3e-10]
field_values = [md.results.HydrologySolution.SedimentHead, md.results.HydrologySolution.SedimentHeadResidual]
