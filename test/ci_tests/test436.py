#Test Name: SquareSheetShelfSteaEnthalpyRheologiesHO
import pyissm
import numpy as np

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/SquareShelf.exp', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetShelf.py')
md = md.extrude(3, 2.)
md = pyissm.model.param.set_flow_equation(md, HO = 'all')
md.cluster.np = 3
md.timestepping.time_step = 0.
md.thermal.isenthalpy = 1
md.initialization.waterfraction = np.zeros((md.mesh.numberofvertices, ))
md.initialization.watercolumn = np.zeros((md.mesh.numberofvertices, ))

# Execute model
field_names = []
field_tolerances = []
field_values = []
for i in ['LliboutryDuval', 'CuffeyTemperate']:
    print(' ')
    print(' ====== Testing rheology law: ' + i + ' = ')

    md.materials.rheology_law = i
    md = pyissm.model.execute.solve(md, 'Steadystate')
    field_names += ['Vx' + i, 'Vy' + i, 'Vz' + i, 'Vel' + i, 'Pressure' + i,
                    'Temperature' + i, 'Waterfraction' + i, 'Enthalpy' + i]
    field_tolerances += [2e-09, 1e-09, 1e-09, 1e-09, 1e-13, 2e-10, 6e-10, 1e-9]
    field_values += [md.results.SteadystateSolution.Vx,
                     md.results.SteadystateSolution.Vy,
                     md.results.SteadystateSolution.Vz,
                     md.results.SteadystateSolution.Vel,
                     md.results.SteadystateSolution.Pressure,
                     md.results.SteadystateSolution.Temperature,
                     md.results.SteadystateSolution.Waterfraction,
                     md.results.SteadystateSolution.Enthalpy]
