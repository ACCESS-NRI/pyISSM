import pyissm


#test reverse scalar vs forward vectorial drivers in ADOLC, using the test3009 setup, equivalent to test109 setup.
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 100000.)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 1

md.autodiff.isautodiff = True
md.toolkits.DefaultAnalysis = pyissm.tools.config.issm_gsl_solver()


#first run scalar reverse mode:
indep = pyissm.model.classes.independent()
indep.name = 'md.geometry.thickness'
indep.type = 'vertex'
indep.nods = md.mesh.numberofvertices
indep.fos_reverse_index = 1
md.autodiff.independents = [indep]

dep = pyissm.model.classes.dependent()
dep.name = 'MaxVel'
dep.type = 'scalar'
dep.fos_reverse_index = 1
md.autodiff.dependents = [dep]
md.autodiff.driver = 'fos_reverse'

md = pyissm.model.execute.solve(md, 'Transient')
#recover jacobian:
jac_reverse = md.results.TransientSolution[0].AutodiffJacobian

#Fields and tolerances to track changes
field_names = ['Jac Reverse']
field_tolerances = [1e-13]
field_values = [jac_reverse]
