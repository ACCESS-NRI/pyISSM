#Test Name: SquareShelfConstrainedBedSlop2d
import pyissm

md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3
md = pyissm.model.execute.solve(md, 'BedSlope')

#Fields and tolerances to track changes
field_names = ['BedSlopeX', 'BedSlopeY']
field_tolerances = [1e-13, 1e-13]
field_values = [md.results.BedSlopeSolution.BedSlopeX,
                md.results.BedSlopeSolution.BedSlopeY]
