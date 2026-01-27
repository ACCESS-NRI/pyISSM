#Test Name: SquareShelfConstrainedSurfSlope3d
import pyissm

md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
md = md.extrude(5, 1.)
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3
md = pyissm.model.execute.solve(md, 'SurfaceSlope')

#Fields and tolerances to track changes
field_names = ['SurfaceSlopeX', 'SurfaceSlopeY']
field_tolerances = [1e-13, 1e-13]
field_values = [md.results.SurfaceSlopeSolution.SurfaceSlopeX,
                md.results.SurfaceSlopeSolution.SurfaceSlopeY]
