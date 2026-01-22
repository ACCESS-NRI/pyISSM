#Test Name: SquareShelfConstrainedBalThic2d
import pyissm

md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
#Add boundary conditions on thickness on the border
pos = np.nonzero(md.mesh.vertexonboundary)
md.balancethickness.spcthickness[pos] = md.geometry.thickness[pos]
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3
md = pyissm.model.execute.solve(md, 'Balancethickness')

#Fields and tolerances to track changes
field_names = ['Thickness']
field_tolerances = [1e-13]
field_values = [md.results.BalancethicknessSolution.Thickness]
