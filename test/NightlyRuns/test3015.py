#Test Name: SquareShelfConstrainedMasstransp2dAdolcForwardDifference
import copy
import pyissm

#This test runs test3015 with autodiff on, and checks that
#the value of the scalar forward difference match a step-wise differential

#First configure
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 50000.)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 1
md.masstransport.requested_outputs = ['IceVolume']
md.toolkits.DefaultAnalysis = pyissm.tools.config.issm_gsl_solver()

#setup autodiff parameters
index = 1  #this is the scalar component we are checking against
md.autodiff.independents = pyissm.model.classes.independent()
md.autodiff.independents.name = md.geometry.thickness
md.autodiff.independents.type = 'vertex'
md.autodiff.independents.nods = md.mesh.numberofvertices
md.autodiff.independents.fos_forward_index = index

md.autodiff.dependents = pyissm.model.classes.dependent()
md.autodiff.dependents.name = 'IceVolume'
md.autodiff.dependents.type = 'scalar'
md.autodiff.driver = 'fos_forward'

#PYTHON: indices start at 0, make sure to offset index
index = index - 1

#parameters for the step-wise derivative
delta = 0.001
h1 = md.geometry.thickness[index]
h0 = h1 * (1. - delta)
h2 = h1 * (1. + delta)
deltaH = (h2 - h0)

#save model
md2 = copy.deepcopy(md)

#evaluate derivative by forward and backward stepping
#forward
md = copy.deepcopy(md2)
md.autodiff.isautodiff = False
md.geometry.thickness[index] = h0
md.geometry.base = -md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness
md.geometry.surface = md.geometry.base + md.geometry.thickness
md = pyissm.model.bc.set_ice_shelf_bc(md)

md = pyissm.model.execute.solve(md, 'Masstransport')
V0 = md.results.MasstransportSolution.IceVolume

#backward
md = copy.deepcopy(md2)
md.autodiff.isautodiff = False
md.geometry.thickness[index] = h2
md.geometry.base = -md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness
md.geometry.surface = md.geometry.base + md.geometry.thickness
md = pyissm.model.classes.bc.set_ice_shelf_bc(md)

md = pyissm.model.execute.solve(md, 'Masstransport')
V2 = md.results.MasstransportSolution.IceVolume

#compute resulting derivative
dVdh_an = (V2 - V0) / deltaH

#evaluate derivative using ADOLC
md = md2
md.autodiff.isautodiff = True
md.geometry.thickness[index] = h1
md.geometry.base = -md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness
md.geometry.surface = md.geometry.base + md.geometry.thickness
md = pyissm.model.bc.set_ice_shelf_bc(md)

md = pyissm.model.execute.solve(md, 'Masstransport')
#retrieve directly
dVdh_ad = md.results.MasstransportSolution.AutodiffJacobian

print("dV / dh: analytical:  %16.16g\n       using adolc:  %16.16g\n" % (dVdh_an, dVdh_ad))

#Fields and tolerances to track changes
field_names = ['dV/dh']
field_tolerances = [1e-8]
field_values = [dVdh_ad]
