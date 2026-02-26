#Test Name: SquareSheetShelfGroundingLine3dSoft
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/SquareShelf.exp', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetShelf.py')
md.cluster.np = 3
md.initialization.vx[:] = 0.
md.initialization.vy[:] = 0.
md.geometry.base = -700. - abs(md.mesh.y - 500000.) / 1000.
md.geometry.bed = -700. - abs(md.mesh.y - 500000.) / 1000.
md.geometry.thickness[:] = 1300
md.geometry.surface = md.geometry.base + md.geometry.thickness
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md = md.extrude(3, 1.)

md.smb.mass_balance[:] = -150
md.transient.isstressbalance = False
md.transient.isgroundingline = True
md.groundingline.migration = 'SoftMigration'
md.groundingline.melt_interpolation = 'SubelementMelt1'

# Execute model
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Bed1', 'Surface1', 'Thickness1', 'Floatingice1',
               'Bed2', 'Surface2', 'Thickness2', 'Floatingice2',
               'Bed3', 'Surface3', 'Thickness3', 'Floatingice3']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13,
                    5e-12, 4e-12, 7e-12, 2e-11,
                    1e-10, 5e-11, 1e-10, 7e-10]
field_values = [md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].MaskOceanLevelset,
                md.results.TransientSolution[1].Base,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].MaskOceanLevelset,
                md.results.TransientSolution[2].Base,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].MaskOceanLevelset]
