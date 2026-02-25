#Test Name: SquareShelfConstrainedStressFSAdolc
import pyissm


md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 180000.)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
md = md.extrude(3, 2.)
md = pyissm.model.param.set_flow_equation(md, FS = 'all')
md.cluster.np = 1
md.autodiff.isautodiff = True
md.toolkits.DefaultAnalysis = pyissm.tools.config.issm_gsl_solver()
md = pyissm.model.execute.solve(md, 'Stressbalance')

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure']
field_tolerances = [1e-08, 1e-08, 3e-07, 1e-08, 1e-08]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vz,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure]
