#Test Name: SquareSheetShelfGroundingLine2dAggressive
import pyissm
import numpy as np

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/SquareShelf.exp', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetShelf.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3
md.initialization.vx[:] = 0.
md.initialization.vy[:] = 0.
md.geometry.base = -700. - abs(md.mesh.y - 500000.) / 1000.
md.geometry.bed = -700. - abs(md.mesh.y - 500000.) / 1000.
md.geometry.thickness[:] = 1000.
md.geometry.surface = md.geometry.base + md.geometry.thickness
md.smb.mass_balance[:] = 100.
md.transient.isstressbalance = False
md.transient.isgroundingline = True
md.groundingline.migration = 'AggressiveMigration'
md.transient.requested_outputs = ['IceVolume', 'IceVolumeAboveFloatation', 'IceVolumeAboveFloatationScaled', 'GroundedArea', 'FloatingArea', 'GroundedAreaScaled', 'FloatingAreaScaled']
md.mesh.scale_factor = 1.1 * np.ones((md.mesh.numberofvertices))

# Execute model
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Bed1', 'Surface1', 'Thickness1', 'Floatingice1', 'IceVolume1', 'IceVolumeAboveFloatation1', 'IceVolumeAboveFloatationScaled1', 'GroundedArea1', 'GroundedAreaScaled1', 'FloatingArea1', 'FloatingAreaScaled1',
               'Bed2', 'Surface2', 'Thickness2', 'Floatingice2', 'IceVolume2', 'IceVolumeAboveFloatation2', 'IceVolumeAboveFloatationScaled2', 'GroundedAred2', 'GroundedAreaScaled2', 'FloatingArea2', 'FloatingAreaScaled2',
               'Bed3', 'Surface3', 'Thickness3', 'Floatingice3', 'IceVolume3', 'IceVolumeAboveFloatation3', 'IceVolumeAboveFloatationScaled3', 'GroundedArea3', 'GroundedAreaScaled3', 'FloatingArea3', 'FloatingAreaScaled3']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].MaskOceanLevelset,
                md.results.TransientSolution[0].IceVolume,
                md.results.TransientSolution[0].IceVolumeAboveFloatation,
                md.results.TransientSolution[0].IceVolumeAboveFloatationScaled,
                md.results.TransientSolution[0].GroundedArea,
                md.results.TransientSolution[0].GroundedAreaScaled,
                md.results.TransientSolution[0].FloatingArea,
                md.results.TransientSolution[0].FloatingAreaScaled,
                md.results.TransientSolution[1].Base,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].MaskOceanLevelset,
                md.results.TransientSolution[1].IceVolume,
                md.results.TransientSolution[1].IceVolumeAboveFloatation,
                md.results.TransientSolution[1].IceVolumeAboveFloatationScaled,
                md.results.TransientSolution[1].GroundedArea,
                md.results.TransientSolution[1].GroundedAreaScaled,
                md.results.TransientSolution[1].FloatingArea,
                md.results.TransientSolution[1].FloatingAreaScaled,
                md.results.TransientSolution[2].Base,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].MaskOceanLevelset,
                md.results.TransientSolution[2].IceVolume,
                md.results.TransientSolution[2].IceVolumeAboveFloatation,
                md.results.TransientSolution[2].IceVolumeAboveFloatationScaled,
                md.results.TransientSolution[2].GroundedArea,
                md.results.TransientSolution[2].GroundedAreaScaled,
                md.results.TransientSolution[2].FloatingArea,
                md.results.TransientSolution[2].FloatingAreaScaled]
