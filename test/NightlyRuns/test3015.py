#Test Name: SquareShelfConstrainedMasstransp2dAdolcForwardDifference
import copy
import pyissm

#This test runs test3015 with autodiff on, and checks that
#the value of the scalar forward difference match a step-wise differential

# First configure
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 50000.)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 1
md.masstransport.requested_outputs = ['IceVolume']
md.toolkits.DefaultAnalysis = pyissm.tools.config.issm_gsl_solver()

# Setup autodiff parameters
index = 1  #this is the scalar component we are checking against
indep = pyissm.model.classes.independent()
indep.name = 'md.geometry.thickness'
indep.type = 'vertex'
indep.nods = md.mesh.numberofvertices
indep.fos_forward_index = index
md.autodiff.independents = [indep]

dep = pyissm.model.classes.dependent()
dep.name = 'IceVolume'
dep.type = 'scalar'
md.autodiff.dependents = [dep]

md.autodiff.driver = 'fos_forward'

# PYTHON: indices start at 0, make sure to offset index
index = index - 1

# Parameters for the step-wise derivative
delta = 0.001
h1 = md.geometry.thickness[index]
h0 = h1 * (1. - delta)
h2 = h1 * (1. + delta)
deltaH = (h2 - h0)

# Save model
md2 = copy.deepcopy(md)

# Evaluate derivative by forward and backward stepping
# forward
md = copy.deepcopy(md2)
md.autodiff.isautodiff = False
md.geometry.thickness[index] = h0
md.geometry.base = -md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness
md.geometry.surface = md.geometry.base + md.geometry.thickness
md = pyissm.model.bc.set_ice_shelf_bc(md)

md = pyissm.model.execute.solve(md, 'Masstransport')
V0 = md.results.MasstransportSolution.IceVolume

# Backward
md = copy.deepcopy(md2)
md.autodiff.isautodiff = False
md.geometry.thickness[index] = h2
md.geometry.base = -md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness
md.geometry.surface = md.geometry.base + md.geometry.thickness
md = pyissm.model.bc.set_ice_shelf_bc(md)

md = pyissm.model.execute.solve(md, 'Masstransport')
V2 = md.results.MasstransportSolution.IceVolume

# Compute resulting derivative
dVdh_an = ((V2 - V0) / deltaH)[0]

# Evaluate derivative using ADOLC
md = md2
md.autodiff.isautodiff = True
md.geometry.thickness[index] = h1
md.geometry.base = -md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness
md.geometry.surface = md.geometry.base + md.geometry.thickness
md = pyissm.model.bc.set_ice_shelf_bc(md)

md = pyissm.model.execute.solve(md, 'Masstransport')

# Retrieve directly
dVdh_ad = md.results.MasstransportSolution.AutodiffJacobian[0][0]

print("dV / dh: analytical:  %16.16g\n       using adolc:  %16.16g\n" % (dVdh_an, dVdh_ad))

#Fields and tolerances to track changes
field_names = ['dV/dh']
field_tolerances = [1e-8]
field_values = [dVdh_ad]