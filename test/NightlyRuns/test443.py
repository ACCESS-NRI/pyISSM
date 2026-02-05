#Test Name: SquareSheetShelfTranSSA3dAggressiveRegionalOutput
import pyissm
import numpy as np
import copy

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/SquareShelf.exp', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetShelf.py')
md.geometry.bed = copy.deepcopy(md.geometry.base)
pos = np.nonzero(md.mask.ocean_levelset < 0.)
md.geometry.bed[pos] = md.geometry.bed[pos] - 10
md.friction = pyissm.model.classes.friction.coulomb()
md.friction.coefficient = 20 * np.ones(md.mesh.numberofvertices)
md.friction.p = 1 * np.ones(md.mesh.numberofelements)
md.friction.q = 1 * np.ones(md.mesh.numberofelements)
md.friction.coefficientcoulomb = 0.02 * np.ones(md.mesh.numberofvertices)
md.transient.isthermal = False
md.transient.isgroundingline = True
md.groundingline.migration = 'AggressiveMigration'
md.mesh.scale_factor = 1.1 * np.ones((md.mesh.numberofvertices))

md.settings.output_frequency = 2
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3

regionalmask = np.zeros((md.mesh.numberofvertices))
inflag = pyissm.tools.wrappers.ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y, '../assets/Exp/SquareHalfRight.exp', 'node', 1)
regionalmask[np.nonzero(inflag)[0]] = 1.
md.transient.requested_outputs = ['default', 'GroundedArea1', 'FloatingArea1', 'TotalFloatingBmb1', 'TotalGroundedBmb1', 'TotalSmb1', 'IceMass1', 'IceVolume1', 'IceVolumeAboveFloatation1', 'IceVolumeAboveFloatation2', 'IceVolumeAboveFloatation', 'IceMassScaled1', 'IceVolumeScaled1', 'IceVolumeAboveFloatationScaled1', 'IceVolumeAboveFloatationScaled2']

md.outputdefinition.definitions = [pyissm.model.classes.regionaloutput(name = 'GroundedArea1', outputnamestring = 'GroundedArea', mask = regionalmask, definitionstring = 'Outputdefinition1'),
                                   pyissm.model.classes.regionaloutput(name = 'FloatingArea1', outputnamestring = 'FloatingArea', mask = regionalmask, definitionstring = 'Outputdefinition2'),
                                   pyissm.model.classes.regionaloutput(name = 'TotalFloatingBmb1', outputnamestring = 'TotalFloatingBmb', mask = regionalmask, definitionstring = 'Outputdefinition3'),
                                   pyissm.model.classes.regionaloutput(name = 'TotalGroundedBmb1', outputnamestring = 'TotalGroundedBmb', mask = regionalmask, definitionstring = 'Outputdefinition4'),
                                   pyissm.model.classes.regionaloutput(name = 'IceMass1', outputnamestring = 'IceMass', mask = regionalmask, definitionstring = 'Outputdefinition5'),
                                   pyissm.model.classes.regionaloutput(name = 'IceVolume1', outputnamestring = 'IceVolume', mask = regionalmask, definitionstring = 'Outputdefinition6'),
                                   pyissm.model.classes.regionaloutput(name = 'IceVolumeAboveFloatation1', outputnamestring = 'IceVolumeAboveFloatation', mask = regionalmask, definitionstring = 'Outputdefinition7'),
                                   pyissm.model.classes.regionaloutput(name = 'TotalSmb1', outputnamestring = 'TotalSmb', mask = regionalmask, definitionstring = 'Outputdefinition8'),
                                   pyissm.model.classes.regionaloutput(name = 'IceVolumeAboveFloatation2', outputnamestring = 'IceVolumeAboveFloatation', maskexpstring = '../assets/Exp/SquareHalfRight.exp', definitionstring = 'Outputdefinition9'),
                                   pyissm.model.classes.regionaloutput(name = 'IceMassScaled1', outputnamestring = 'IceMassScaled', mask = regionalmask, definitionstring = 'Outputdefinition10'),
                                   pyissm.model.classes.regionaloutput(name = 'IceVolumeScaled1', outputnamestring = 'IceVolumeScaled', mask = regionalmask,  definitionstring = 'Outputdefinition11'),
                                   pyissm.model.classes.regionaloutput(name = 'IceVolumeAboveFloatationScaled1', outputnamestring = 'IceVolumeAboveFloatationScaled', mask = regionalmask, definitionstring = 'Outputdefinition12'),
                                   pyissm.model.classes.regionaloutput(name = 'IceVolumeAboveFloatationScaled2', outputnamestring = 'IceVolumeAboveFloatationScaled', maskexpstring = '../assets/Exp/SquareHalfRight.exp', definitionstring = 'Outputdefinition13')]

md = md.extrude(3, 1.)
md2 = copy.deepcopy(md)
md2.collapse()
md.groundingline.intrusion_distance = np.zeros(md.mesh.numberofvertices, )

# Execute model
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['IceMass1', 'IceVolume1', 'IceVolumeAboveFloatation1', 'IceVolumeAboveFloatation21', 'Thickness1', 'GroundedArea1', 'FloatingArea1', 'TotalFloatingBmb1', 'TotalGroundedBmb1', 'TotalSmb1', 'IceMassScaled1', 'IceVolumeScaled1', 'IceVolumeAboveFloatationScaled1', 'IceVolumeAboveFloatationScaled21', 'IceMass3', 'IceVolume3', 'IceVolumeAboveFloatation3', 'IceVolumeAboveFloatation23', 'Thickness3', 'GroundedArea3', 'FloatingArea3', 'TotalFloatingBmb3', 'TotalGroundedBmb3', 'TotalSmb3', 'IceMassScaled3', 'IceVolumeScaled3', 'IceVolumeAboveFloatationScaled3', 'IceVolumeAboveFloatationScaled23', 'ExtrudedMask', 'CollapsedMask']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md.results.TransientSolution[0].IceMass1,
                md.results.TransientSolution[0].IceVolume1,
                md.results.TransientSolution[0].IceVolumeAboveFloatation1,
                md.results.TransientSolution[0].IceVolumeAboveFloatation2,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].GroundedArea1,
                md.results.TransientSolution[0].FloatingArea1,
                md.results.TransientSolution[0].TotalFloatingBmb1,
                md.results.TransientSolution[0].TotalGroundedBmb1,
                md.results.TransientSolution[0].TotalSmb1,
                md.results.TransientSolution[0].IceMassScaled1,
                md.results.TransientSolution[0].IceVolumeScaled1,
                md.results.TransientSolution[0].IceVolumeAboveFloatationScaled1,
                md.results.TransientSolution[0].IceVolumeAboveFloatationScaled2,
                md.results.TransientSolution[2].IceMass1,
                md.results.TransientSolution[2].IceVolume1,
                md.results.TransientSolution[2].IceVolumeAboveFloatation1,
                md.results.TransientSolution[2].IceVolumeAboveFloatation2,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].GroundedArea1,
                md.results.TransientSolution[2].FloatingArea1,
                md.results.TransientSolution[2].TotalFloatingBmb1,
                md.results.TransientSolution[2].TotalGroundedBmb1,
                md.results.TransientSolution[2].TotalSmb1,
                md.results.TransientSolution[2].IceMassScaled1,
                md.results.TransientSolution[2].IceVolumeScaled1,
                md.results.TransientSolution[2].IceVolumeAboveFloatationScaled1,
                md.results.TransientSolution[2].IceVolumeAboveFloatationScaled2,
                len(md.outputdefinition.definitions[0].mask),
                len(md2.outputdefinition.definitions[0].mask)]
