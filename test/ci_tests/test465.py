#Test Name: SquareSheetShelfAmrBamgAll
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/SquareShelf.exp', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetShelf.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3
md.transient.isstressbalance = 1
md.transient.ismasstransport = 1
md.transient.issmb = 0
md.transient.isthermal = 0
md.transient.isgroundingline = 0
# amr bamg settings, field, grounding line and ice front
md.amr.hmin = 20000
md.amr.hmax = 100000
md.amr.fieldname = 'Vel'
md.amr.keepmetric = 0
md.amr.gradation = 1.2
md.amr.groundingline_resolution = 10000
md.amr.groundingline_distance = 100000
md.amr.icefront_resolution = 15000
md.amr.icefront_distance = 100000
md.amr.thicknesserror_resolution = 1000
md.amr.thicknesserror_threshold = 0
md.amr.deviatoricerror_resolution = 1000
md.amr.deviatoricerror_threshold = 0
md.transient.amr_frequency = 1
md.timestepping.start_time = 0
md.timestepping.final_time = 3
md.timestepping.time_step = 1

# Execute model
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vel', 'Pressure']
field_tolerances = [1e-8, 1e-8, 1e-8, 1e-8]
field_values = [md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Vel,
                md.results.TransientSolution[2].Pressure]
