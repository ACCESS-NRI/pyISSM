#Test Name: SquareShelfConstrainedStressSSA2dAdolcMumps
import pyissm


md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 50000.)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3
md.stressbalance.requested_outputs = ['default', 'DeviatoricStressxx', 'DeviatoricStressyy', 'DeviatoricStressxy']

md.autodiff.isautodiff = True
md.toolkits.DefaultAnalysis = pyissm.tools.config.issm_mumps_solver()
md = pyissm.model.execute.solve(md, 'Stressbalance')

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vel', 'Pressure',
               'DeviatoricStressxx', 'DeviatoricStressyy', 'DeviatoricStressxy']
field_tolerances = [1e-12, 1e-12, 1e-12, 1e-12,
                    1e-12, 1e-12, 1e-12]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure,
                md.results.StressbalanceSolution.DeviatoricStressxx,
                md.results.StressbalanceSolution.DeviatoricStressyy,
                md.results.StressbalanceSolution.DeviatoricStressxy]
