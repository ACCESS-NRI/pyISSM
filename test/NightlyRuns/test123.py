#Test Name: SquareShelfConstrainedTranMisfitSurface
import pyissm
import numpy as np

md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 180000)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
md = pyissm.model.param.set_flow_equation(md, SSA =  'all')
md.cluster.np = 3

fake_surface = np.vstack((np.append(np.array(md.geometry.surface) + 100, 1.1),
                          np.append(np.array(md.geometry.surface) + 200, 2.1),
                          np.append(np.array(md.geometry.surface) + 300, 2.5))).T

md.misfit = pyissm.model.classes.misfit()
md.misfit.name = 'SurfaceMisfit'
md.misfit.definitionstring = 'Outputdefinition1'
md.misfit.model_string = 'Surface'
md.misfit.observation = fake_surface
md.misfit.observation_string = 'SurfaceObservation'
md.misfit.timeinterpolation = 'nearestneighbor'
md.misfit.weights = np.ones((md.mesh.numberofvertices, 1))
md.misfit.weights_string = 'WeightsSurfaceObservation'

md.transient.requested_outputs = ['default', 'SurfaceMisfit']
md.outputdefinition.definitions = [md.misfit]

md = pyissm.model.execute.solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['SurfaceMisfitFirstStep', 'SurfaceMisfitSecondStep', 'SurfaceMisfitThirdStep']
field_tolerances = [1e-13, 1e-13, 1e-13]
field_values = [md.results.TransientSolution[0].SurfaceMisfit,
                md.results.TransientSolution[1].SurfaceMisfit,
                md.results.TransientSolution[2].SurfaceMisfit]
