#Test Name: SquareShelfConstrainedTherTranAdolcReverseVsForward
import numpy as np
import pyissm

# test reverse scalar vs forward vectorial drivers in ADOLC, using the test3009 setup, equivalent to test109 setup.
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 100000.)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 1
md.toolkits.DefaultAnalysis = pyissm.tools.config.issm_gsl_solver()

md.autodiff.isautodiff = True
md.verbose.autodiff = True

# first run scalar reverse mode:
indep = pyissm.model.classes.independent()
indep.name = 'md.geometry.thickness'
indep.type = 'vertex'
indep.nods = md.mesh.numberofvertices
indep.fos_reverse_index = 1
md.autodiff.independents = [indep]

dep = pyissm.model.classes.dependent()
dep.name = 'MaxVel'
dep.type = 'scalar'
dep.nods = md.mesh.numberofvertices
dep.fos_reverse_index = 1
md.autodiff.dependents = [dep]

md.autodiff.driver = 'fos_reverse'

md = pyissm.model.execute.solve(md, 'Transient')

# recover jacobian:
jac_reverse = md.results.TransientSolution[0].AutodiffJacobian

# now run vectorial forward mode
indep = pyissm.model.classes.independent()
indep.name = 'md.geometry.thickness'
indep.type = 'vertex'
indep.nods = md.mesh.numberofvertices
# "not empty" so the MATLAB branch triggers:
indep.fos_forward_index = 1

# scalar “how many directions” (not the explicit list 1..N):
indep.fov_forward_indices = int(md.mesh.numberofvertices)
md.autodiff.independents = [indep]

dep = pyissm.model.classes.dependent()
dep.name = 'MaxVel'
dep.type = 'scalar'
dep.nods = md.mesh.numberofvertices
md.autodiff.dependents = [dep]

md.autodiff.driver = 'fov_forward'

# def debug_isnan_compat(obj, label="obj"):
#     bad = []
#     for k, v in obj.__dict__.items():
#         try:
#             a = np.asarray(v)
#             if a.dtype == object:
#                 bad.append((k, type(v), "dtype=object"))
#             else:
#                 # try isnan if it looks numeric-ish
#                 np.isnan(a)
#         except Exception as e:
#             bad.append((k, type(v), str(e)))
#     if bad:
#         print(f"\n[{label}] potential isnan offenders:")
#         for item in bad:
#             print(" ", item)
#     else:
#         print(f"\n[{label}] looks isnan-safe")


# debug_isnan_compat(md.autodiff, "md.autodiff")
# debug_isnan_compat(md.autodiff.independents, "md.autodiff.independents")
# debug_isnan_compat(md.autodiff.dependents, "md.autodiff.dependents")
md = pyissm.model.execute.solve(md, 'Transient')

# recover jacobian:
jac_forward = md.results.TransientSolution[0].AutodiffJacobian

# Fields and tolerances to track changes
field_names = ['Jac Forward', 'Jac Reverse', 'Jac Forward - Reverse']
field_tolerances = [1e-8, 1e-8, 5e-6]
field_values = [jac_forward, jac_reverse, jac_forward - jac_reverse]