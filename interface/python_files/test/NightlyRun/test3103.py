#Test Name: SquareShelfConstrainedStressHOAdolcMumps
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from issmmumpssolver import issmmumpssolver


md = triangle(model(), '../Exp/Square.exp', 180000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelfConstrained.py')
md.extrude(3, 2.)
md = setflowequation(md, 'HO', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md.stressbalance.requested_outputs = ['default', 'StressTensorxx', 'StressTensoryy', 'StressTensorzz', 'StressTensorxy', 'StressTensorxz', 'StressTensoryz']

md.autodiff.isautodiff = True
md.toolkits.DefaultAnalysis = issmmumpssolver()
md = solve(md, 'Stressbalance')

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure',
               'StressTensorxx', 'StressTensoryy', 'StressTensorzz', 'StressTensorxy', 'StressTensorxz', 'StressTensoryz']
field_tolerances = [1e-09, 1e-09, 1e-09, 1e-09, 1e-09,
                    1e-09, 1e-09, 1e-09, 1e-09, 1e-09, 2e-09]
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
