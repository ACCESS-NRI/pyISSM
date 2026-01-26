#Test Name: SquareSheetHydrologyGlaDSSheetInPhi
import pyissm
import numpy as np

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 50000.)
md.mesh.x = md.mesh.x / 100
md.mesh.y = md.mesh.y / 100
md = pyissm.model.param.set_mask(md, None, None)
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.stressbalance.maxiter = 10
md.cluster.np = 2

# Some constants
md.constants.g = 9.8
md.materials.rho_ice = 910

# Geometry
md.geometry.surface = -0.02 * md.mesh.x + 320
md.geometry.bed = np.zeros((md.mesh.numberofvertices))
md.geometry.base = md.geometry.bed
md.geometry.thickness = md.geometry.surface - md.geometry.bed

# Define initial conditions
md.initialization.vx = 1.0e-6 * md.constants.yts * np.ones((md.mesh.numberofvertices))
md.initialization.vy = np.zeros((md.mesh.numberofvertices))
md.initialization.temperature = (273. - 20.) * np.ones((md.mesh.numberofvertices))
md.initialization.watercolumn = 0.03 * np.ones((md.mesh.numberofvertices))
md.initialization.hydraulic_potential = md.materials.rho_ice * md.constants.g * md.geometry.thickness

# Materials
md.materials.rheology_B = (5e-25)**(-1./3.) * np.ones((md.mesh.numberofvertices))
md.materials.rheology_n = 3. * np.ones((md.mesh.numberofelements))

# Friction
md.friction.coefficient = np.zeros((md.mesh.numberofvertices))
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))
#md.friction.coupling = 0

# Boundary conditions:
md = pyissm.model.bc.set_ice_sheet_bc(md)

md.inversion.iscontrol = 0
md.transient = pyissm.model.classes.transient.deactivate_all(md.transient)
md.transient.ishydrology = 1

# Set numerical conditions
md.timestepping.time_step = 0.1 / 365
md.timestepping.final_time = 0.4 / 365

# Change hydrology class to Glads model
md.hydrology = pyissm.model.classes.hydrology.glads()
md.hydrology.ischannels = 1
md.hydrology.isincludesheetthickness = 1
md.hydrology.englacial_void_ratio = 1.e-5
md.hydrology.moulin_input = np.zeros((md.mesh.numberofvertices))
md.hydrology.neumannflux = np.zeros((md.mesh.numberofelements))
md.hydrology.bump_height = 1.e-1 * np.ones((md.mesh.numberofvertices))
md.hydrology.sheet_conductivity = 1.e-3 * np.ones((md.mesh.numberofvertices))
md.hydrology.channel_conductivity = 5.e-2 * np.ones((md.mesh.numberofvertices))
md.hydrology.rheology_B_base = 8.378836055370960e+07 * np.ones((md.mesh.numberofvertices)) 

# BCs for hydrology
pos = np.where(np.logical_and(md.mesh.x == 100, md.mesh.vertexonboundary))
md.hydrology.spcphi = np.nan * np.ones((md.mesh.numberofvertices))
md.hydrology.spcphi[pos] = md.materials.rho_ice * md.constants.g * md.geometry.thickness[pos]

# Execute model
md.miscellaneous.name = 'testChannels'
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['HydrologySheetThickness1', 'HydraulicPotential1', 'ChannelArea1',
               'HydrologySheetThickness2', 'HydraulicPotential2', 'ChannelArea2',
               'HydrologySheetThickness3', 'HydraulicPotential3', 'ChannelArea3',
               'HydrologySheetThickness4', 'HydraulicPotential4', 'ChannelArea4']
field_tolerances = [1e-14, 7e-14, 3e-12,
                    1e-14, 7e-14, 3e-12,
                    1e-14, 7e-14, 3e-12,
                    1e-14, 8e-14, 3e-12]
field_values = [md.results.TransientSolution[0].HydrologySheetThickness,
                md.results.TransientSolution[0].HydraulicPotential,
                md.results.TransientSolution[0].ChannelArea,
                md.results.TransientSolution[1].HydrologySheetThickness,
                md.results.TransientSolution[1].HydraulicPotential,
                md.results.TransientSolution[1].ChannelArea,
                md.results.TransientSolution[2].HydrologySheetThickness,
                md.results.TransientSolution[2].HydraulicPotential,
                md.results.TransientSolution[2].ChannelArea,
                md.results.TransientSolution[3].HydrologySheetThickness,
                md.results.TransientSolution[3].HydraulicPotential,
                md.results.TransientSolution[3].ChannelArea]
