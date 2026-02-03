#Test Name: ThermalConduction
import numpy as np
import pyissm


"""
This file can be run to check that the conduction is correctly modeled.
There is no velocity (no advection) the only thermal boundary conditions are an imposed temperature
at the lower and upper surface. The result must be a linear temperature from the upper to the lower
surface. if it is not the case, something is thermal modeling has been changed...
"""

printingflag = False

md = pyissm.model.Model()
md = pyissm.model.mesh.triangle(md, '../assets/Exp/Square.exp', 100000.)
md = pyissm.model.param.set_mask(md, 'all', '')
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareThermal.py')
md = md.extrude(11, 2.)
md = pyissm.model.param.set_flow_equation(md, HO = 'all')


pos1 = np.where(np.isnan(md.mesh.lowerelements))[0]
md.thermal.spctemperature[md.mesh.elements[pos1, 0:3] - 1] = 10.
pos2 = np.where(np.isnan(md.mesh.upperelements))[0]
md.thermal.spctemperature[md.mesh.elements[pos2, 3:6] - 1] = 0.
md.initialization.pressure = np.zeros((md.mesh.numberofvertices), int)

#analytical results
#d2T / dz2 = 0 T(bed)=10 T(surface)=0 = > T = 0 * (z - bed) / thickness + 10 * (surface-z) / thickness
#each layer of the 3d mesh must have a constant value
md.initialization.temperature = 10. * (md.geometry.surface - md.mesh.z) / md.geometry.thickness

#modeled results
md.cluster.np = 2
md = pyissm.model.execute.solve(md, 'Thermal')

#plot results
comp_temp = md.results.ThermalSolution.Temperature
relative = np.abs((comp_temp - md.initialization.temperature) / md.initialization.temperature) * 100.
relative[np.nonzero(comp_temp == md.initialization.temperature)[0]] = 0.
#plotmodel(md, 'data', comp_temp, 'title', 'Modeled temperature [K]', 'data', md.initialization.temperature, 'view', 3, ...
#       'title', 'Analytical temperature [K]', 'view', 3, 'data', comp_temp - md.initialization.temperature, ...
#       'title', 'Absolute error [K]', 'view', 3, 'data', relative, 'title', 'Relative error [%]', 'view', 3, ...
#       'figposition', 'mathieu', 'FontSize  #all', 20)
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('thermalconduction', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 0.7, 'hardcopy', 'off')
#       system(['mv thermalconduction.png ' ISSM_DIR '/website/doc_pdf/validation/Images/Thermal '])

#Fields and tolerances to track changes
field_names = ['ConductionTemperature']
field_tolerances = [1e-13]
field_values = [comp_temp]
