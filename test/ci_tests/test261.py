#Test Name: SquareShelfConstrainedTranEnhanced
import numpy as np
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 180000)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
md = md.extrude(3, 1.)
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3
md.materials = pyissm.model.classes.materials.enhancedice()
md.materials.rheology_B = 3.15e8 * np.ones(md.mesh.numberofvertices, )
md.materials.rheology_n = 3 * np.ones(md.mesh.numberofelements, )
md.materials.rheology_E = np.ones(md.mesh.numberofvertices, )
md.transient.isstressbalance = 1
md.transient.ismasstransport = 0
md.transient.issmb = 1
md.transient.isthermal = 1
md.transient.isgroundingline = 0

# Execute model
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vel', 'Temperature', 'BasalforcingsGroundediceMeltingRate']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Temperature,
                md.results.TransientSolution[0].BasalforcingsGroundediceMeltingRate]
