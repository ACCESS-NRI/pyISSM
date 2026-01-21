#Test Name: SquareShelfConstrainedStressHO
import pyissm

md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 180000)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
md = md.extrude(3, 2.)
md = pyissm.model.param.set_flow_equation(md, 'HO', 'all')
md.cluster.np = 3
md.stressbalance.requested_outputs = ['default', 'StressTensorxx', 'StressTensoryy', 'StressTensorzz', 'StressTensorxy', 'StressTensorxz', 'StressTensoryz']
md = pyissm.model.execute.solve(md, 'Stressbalance')
#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz',
               'Vel', 'Pressure',
               'StressTensorxx', 'StressTensoryy', 'StressTensorzz',
               'StressTensorxy', 'StressTensorxz', 'StressTensoryz']
field_tolerances = [1e-09, 1e-09, 1e-09,
                    1e-09, 1e-09,
                    1e-09, 1e-09, 1e-09,
                    1e-09, 1e-09, 1e-08]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vz,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure,
                md.results.StressbalanceSolution.StressTensorxx,
                md.results.StressbalanceSolution.StressTensoryy,
                md.results.StressbalanceSolution.StressTensorzz,
                md.results.StressbalanceSolution.StressTensorxy,
                md.results.StressbalanceSolution.StressTensorxz,
                md.results.StressbalanceSolution.StressTensoryz]
