#Test Name: SquareSheetConstrainedGia2d
#GIA test, based off of test101. Running default GIA Ivins class.
import numpy as np
from socket import gethostname
import pyissm


md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 100000.)
md = pyissm.model.param.set_mask(md, None, None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetConstrained.py')

# GIA Ivins, 2 layer model
md.solidearth.settings.grdmodel = 2
md.solidearth.settings.isgrd = 1

md.materials = pyissm.model.classes.materials('litho','ice')
md.materials.numlayers = 2;
md.materials.radius = [10, 6271e3, 6371e3]
md.materials.density = [3.34e3, 3.32e3]
md.materials.lame_mu = [1.45e11, 6.7e10]
md.materials.viscosity = [1e21, 0]
md.initialization.sealevel = np.zeros(md.mesh.numberofvertices)
md.solidearth.settings.cross_section_shape = 1 # for square-edged x-section
md.solidearth.settings.grdocean = 0 # do not compute sea level, only deformation
md.solidearth.settings.sealevelloading = 0 # do not compute sea level, only deformation
md.solidearth.requested_outputs = ['Sealevel', 'BedGRD']

# Loading history
md.timestepping.start_time = -2400000  # 4,800 kyr :: EVALUATION TIME
md.timestepping.time_step = 2400000  # 2,400 kyr :: EVALUATION TIME
# To get rid of default final_time, make sure final_time > start_time
md.timestepping.final_time = 2400000  # 2,400 kyr
md.masstransport.spcthickness = np.array([
    np.append(md.geometry.thickness, 0),
    np.append(md.geometry.thickness, 2400000)
    ]).T

# Geometry at 0 initially
md.geometry.thickness = np.zeros(md.mesh.numberofvertices)
md.geometry.surface = np.zeros(md.mesh.numberofvertices)
md.geometry.base = np.zeros(md.mesh.numberofvertices)

# Physics
md.transient.issmb = 0
md.transient.isstressbalance = 0
md.transient.isthermal = 0
md.transient.ismasstransport = 1
md.transient.isslc = 1

#Solve for GIA deflection
md.cluster.np = 3
md.verbose.solver = 0
md = pyissm.model.execute.solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['UGrd']
field_tolerances = [1e-13]
field_values = [md.results.TransientSolution[1].BedGRD]
