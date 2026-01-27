#Test Name: 79NorthPISMhydro
import pyissm
import numpy as np

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/79North.exp', 10000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/79NorthShelf.exp', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/79North.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3

# Hydrology
md.hydrology = pyissm.model.classes.hydrology.pism()
md.hydrology.drainage_rate = 0.001 * np.ones((md.mesh.numberofvertices))
md.hydrology.watercolumn_max = 20 * np.ones((md.mesh.numberofvertices))
md.initialization.pressure = np.zeros((md.mesh.numberofvertices))
md.initialization.watercolumn = np.zeros((md.mesh.numberofvertices))
md.basalforcings.groundedice_melting_rate = np.arange(0, md.mesh.numberofvertices) + 1
md.transient.ishydrology = 1
md.transient.issmb = 0
md.transient.ismasstransport = 0
md.transient.isstressbalance = 0
md.transient.isthermal = 0
md.transient.isgroundingline = 0

# Execute model
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['WaterColumn1', 'WaterColumn2', 'WaterColumn3']
field_tolerances = [1e-12, 1e-12, 1e-12]
field_values = [md.results.TransientSolution[0].Watercolumn,
                md.results.TransientSolution[1].Watercolumn,
                md.results.TransientSolution[2].Watercolumn]
