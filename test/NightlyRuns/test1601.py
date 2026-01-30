#Test Name: SquareShelfSSA2dRotation
import numpy as np
import sys
import pyissm


md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000.)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '.../assets/Par/SquareShelf.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.stressbalance.spcvx[np.where(md.mesh.y > 0.)] = np.nan
md.initialization.vx[:] = 0.
md.initialization.vy[:] = 0.
md.initialization.vel = np.zeros((md.mesh.numberofvertices))

md.cluster.np = 2
md = pyissm.model.execute.solve(md, 'Stressbalance')
vel0 = md.results.StressbalanceSolution.Vel
theta = 30. * np.pi / 180.
x = md.mesh.x
y = md.mesh.y
md.mesh.x = np.cos(theta) * x - np.sin(theta) * y
md.mesh.y = np.sin(theta) * x + np.cos(theta) * y

rotation_array = np.array([np.cos(theta), np.sin(theta), 0])
md.stressbalance.referential[:, 0:3] = (np.tile(rotation_array, (md.mesh.numberofvertices, 1)))
md.stressbalance.referential[:, 3:] = np.tile([0, 0, 1], (md.mesh.numberofvertices, 1))
md = pyissm.model.execute.solve(md, 'Stressbalance')
vel1 = md.results.StressbalanceSolution.Vel
#plotmodel(md, 'data', vel0, 'data', vel1, 'data', vel1 - vel0, 'title', 'Cartesian CS', 'title', 'Rotated CS', 'title', 'difference')
print("Error between Cartesian and rotated CS: {}".format(np.max(np.abs(vel0 - vel1)) / (np.max(np.abs(vel0)) + sys.float_info.epsilon)))

#Now, put CS back to normal except on the side where the spc are applied
pos = np.where(np.logical_or(x == 0., x == 1000000.))[0]
md.stressbalance.referential[:] = np.nan
md.stressbalance.referential[pos, 0:3] = np.tile([np.cos(theta), np.sin(theta), 0], (len(pos), 1))
md.stressbalance.referential[pos, 3:] = np.tile([0, 0, 1], (len(pos), 1))
md = pyissm.model.execute.solve(md, 'Stressbalance')
vel2 = md.results.StressbalanceSolution.Vel

#plotmodel(md, 'data', vel0, 'data', vel2, 'data', vel2 - vel0, 'title', 'Cartesian CS', 'title', 'Rotated CS', 'title', 'difference')
print("Error between Cartesian and rotated CS: {}".format(np.max(np.abs(vel0 - vel2)) / (np.max(np.abs(vel0)) + sys.float_info.epsilon)))

#Fields and tolerances to track changes
field_names = ['vel1', 'vel2']
field_tolerances = [1e-11, 1e-11]
field_values = [vel1, vel2]
