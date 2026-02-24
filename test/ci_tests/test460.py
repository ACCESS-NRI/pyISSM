#Test Name: SquareSheetShelfStressFSEstar
import pyissm
import numpy as np

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 180000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/SquareShelf.exp', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetShelf.py')
md = md.extrude(3, 1.)
md.cluster.np = 3
md.materials = pyissm.model.classes.materials.estar()
md.materials.rheology_B = 3.15e8 * np.ones((md.mesh.numberofvertices, ))
md.materials.rheology_Ec = np.ones((md.mesh.numberofvertices, ))
md.materials.rheology_Es = 3 * np.ones((md.mesh.numberofvertices, ))

# Execute model
field_names = []
field_tolerances = []
field_values = []
for i in ['SSA', 'HO', 'FS']:
    if i == 'SSA':
        md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
    elif i == 'HO':
        md = pyissm.model.param.set_flow_equation(md, HO = 'all')
    elif i == 'FS':
        md = pyissm.model.param.set_flow_equation(md, FS = 'all')
    
    md = pyissm.model.execute.solve(md, 'Stressbalance')
    field_names = field_names + ['Vx' + i, 'Vy' + i, 'Vz' + i, 'Vel' + i, 'Pressure' + i]
    field_tolerances = field_tolerances + [7e-06, 2e-05, 2e-06, 5e-06, 8e-07]
    field_values = field_values + [md.results.StressbalanceSolution.Vx,
                                   md.results.StressbalanceSolution.Vy,
                                   md.results.StressbalanceSolution.Vz,
                                   md.results.StressbalanceSolution.Vel,
                                   md.results.StressbalanceSolution.Pressure]
