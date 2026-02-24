#Test Name: ValleyGlacierLevelsetThermalSSA3d
import pyissm
import numpy as np

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 50000)
md = pyissm.model.param.set_mask(md, None, None)
md = pyissm.model.param.parameterize(md, '../assets/Par/ValleyGlacierShelf.py')
md = md.extrude(3, 2.)
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3

# Thermal model
pos_surf = np.nonzero(md.mesh.vertexonsurface)[0]
md.thermal.spctemperature[pos_surf] = md.initialization.temperature[pos_surf]

# Transient
md.transient.isstressbalance = True
md.transient.ismovingfront = True
md.transient.ismasstransport = True
md.transient.issmb = True
md.transient.isthermal = True
md.transient.isgroundingline = True
md.groundingline.melt_interpolation = 'SubelementMelt1'

# Execute model
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'Pressure1', 'Thickness1', 'Surface1', 'MaskIceLevelset1', 'Temperature1',
               'Vx2', 'Vy2', 'Vel2', 'Pressure2', 'Thickness2', 'Surface2', 'MaskIceLevelset2', 'Temperature2',
               'Vx3', 'Vy3', 'Vel3', 'Pressure3', 'Thickness3', 'Surface3', 'MaskIceLevelset3', 'Temperature3']
field_tolerances = [1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12,
                    1e-9, 1e-9, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-9,
                    1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].MaskIceLevelset,
                md.results.TransientSolution[0].Temperature,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Vel,
                md.results.TransientSolution[1].Pressure,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].MaskIceLevelset,
                md.results.TransientSolution[1].Temperature,
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Vel,
                md.results.TransientSolution[2].Pressure,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].MaskIceLevelset,
                md.results.TransientSolution[2].Temperature]
