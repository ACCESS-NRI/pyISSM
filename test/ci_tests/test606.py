#Test Name: 79NorthBedSlop2d
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/79North.exp', 10000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/79NorthShelf.exp', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/79North.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3

# Execute model
md = pyissm.model.execute.solve(md, 'BedSlope')

# Fields and tolerances to track changes
field_names = ['BedSlopeX', 'BedSlopeY']
field_tolerances = [1e-13, 1e-13]
field_values = [md.results.BedSlopeSolution.BedSlopeX,
                md.results.BedSlopeSolution.BedSlopeY]
