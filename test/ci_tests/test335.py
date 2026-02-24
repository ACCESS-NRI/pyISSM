#Test Name: SquareSheetConstrainedStressMOLHO2d
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 200000.)
md = pyissm.model.param.set_mask(md, None, None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetConstrained.py')
md = pyissm.model.param.set_flow_equation(md, MOLHO = 'all')
md = pyissm.model.bc.set_molho_bc(md)
md = md.extrude(5, 1.)
md.cluster.np = 3
md.stressbalance.requested_outputs = ['default', 'VxSurface', 'VySurface', 'VxShear', 'VyShear', 'VxBase', 'VyBase']

# Execute model
md = pyissm.model.execute.solve(md, 'Stressbalance')

# Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure', 'VxSurface', 'VySurface', 'VxShear', 'VyShear', 'VxBase', 'VyBase']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vz,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure,
                md.results.StressbalanceSolution.VxSurface,
                md.results.StressbalanceSolution.VySurface,
                md.results.StressbalanceSolution.VxShear,
                md.results.StressbalanceSolution.VyShear,
                md.results.StressbalanceSolution.VxBase,
                md.results.StressbalanceSolution.VyBase]
