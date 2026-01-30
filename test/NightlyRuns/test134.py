#Test Name: SquareShelfConstrainedSampling
import numpy as np
import pyissm

md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')


md = md.sampling.setparameters(md, 2e5, 1.)
md.sampling.seed = 1
n = md.mesh.numberofvertices
md.sampling.phi = np.zeros((n, 1))
md.initialization.sample = np.ones((n, 1))

md = pyissm.model.execute.solve(md, 'smp')

#Fields and tolerances to track changes
field_names = ['Sample']
field_tolerances = [1e-13]
field_values = [md.results.SamplingSolution.Sample]
