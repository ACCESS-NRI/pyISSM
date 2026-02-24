#Test Name:79NorthHydrologyArmapw
import pyissm
import numpy as np

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/79North.exp', 6000)
md = pyissm.model.param.set_mask(md, '../assets/Exp/79NorthShelf.exp', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/79North.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3

# Friction
md.friction.coefficient = 30 * np.ones(md.mesh.numberofvertices)
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))

# Basin separation default
idb_df = np.zeros((md.mesh.numberofelements))
iid1 = np.where(md.mesh.y <= -1.08e6)[0]
for ii in range(md.mesh.numberofelements):
    for vertex in range(3):
        if md.mesh.elements[ii][vertex] - 1 in iid1:  # one vertex in basin 1; NOTE: offset because of 1-based vertex indexing
            idb_df[ii] = 1
    if idb_df[ii] == 0: # no vertex was found in basin 1
        for vertex in range(3):
            idb_df[ii] = 2
# Covariance matrix
covGlob = np.array([[1e9,0,0,0],[0,1e9,0,0],[0,0,0.1,0],[0,0,0,0.1]])

# Hydrology scheme
md.hydrology = pyissm.model.classes.hydrology.armapw()
md.hydrology.num_basins = 2
md.hydrology.basin_id = np.copy(idb_df).astype(int)
md.hydrology.monthlyfactors = 1 * np.ones((md.hydrology.num_basins,12))
md.hydrology.monthlyfactors[:,0:3] = 0
md.hydrology.num_params = 2 # number of parameters in the polynomial
md.hydrology.num_breaks = 2 # number of breakpoints
termconst = np.array([[0.5 * 1e6, 0.1 * 1e6, 0.5e6],[0.5 * 1e6, 0.1 * 1e6, 0.5e6]])
termtrend = np.array([[1 * 1e5, 0, 0],[0, 1 * 1e5, 0]])
md.hydrology.polynomialparams = np.transpose(np.stack((termconst, termtrend)), (1, 2, 0))
md.hydrology.datebreaks = np.array([[20, 40], [20, 40]])
md.hydrology.arma_timestep = 1
md.hydrology.ar_order = 1
md.hydrology.ma_order = 1
md.hydrology.arlag_coefs = np.array([[0.98], [0.98]])
md.hydrology.malag_coefs = np.array([[0], [0]])

# SMB
md.smb = pyissm.model.classes.smb.arma()
md.smb.num_basins = 2
md.smb.basin_id = np.copy(idb_df)
md.smb.num_breaks = 0
md.smb.num_params = 1
md.smb.polynomialparams = 0 * np.array([[0.5], [0.2]])
md.smb.ar_order = 1
md.smb.ma_order = 1
md.smb.arlag_coefs = np.array([[0], [0]])
md.smb.malag_coefs = np.array([[0], [0]])
md.smb.arma_timestep = 1.0

# Stochastic forcing
md.stochasticforcing.isstochasticforcing = 1
md.stochasticforcing.fields = ['FrictionWaterPressure','SMBarma']
md.stochasticforcing.defaultdimension = 2
md.stochasticforcing.default_id = idb_df
md.stochasticforcing.covariance = covGlob # global covariance
md.stochasticforcing.stochastictimestep  = 1; # time step of stochastic forcing
md.stochasticforcing.randomflag = 0  # determines true/false randomness

md.transient.issmb = 1
md.transient.ismasstransport = 1
md.transient.isstressbalance = 1
md.transient.isthermal = 0
md.transient.isgroundingline = 0
md.transient.ishydrology = 1

md.transient.requested_outputs = ['default','SmbMassBalance', 'FrictionWaterPressure']
md.timestepping.start_time = 0
md.timestepping.time_step = 1.0/12
md.timestepping.final_time = 2

# Execute model
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Vel1',  'Thickness1', 'SmbMassBalance1',  'FrictionWaterPressure1',
               'Vel12', 'Thickness12','SmbMassBalance12', 'FrictionWaterPressure12',
               'Vel24', 'Thickness24','SmbMassBalance24', 'FrictionWaterPressure24']

field_tolerances = [2e-10, 2e-10, 2e-10, 2e-10,
                    4e-10, 4e-10, 4e-10, 4e-10,
                    8e-10, 8e-10, 8e-10, 8e-10]

field_values = [md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].SmbMassBalance,
                md.results.TransientSolution[0].FrictionWaterPressure,
                md.results.TransientSolution[11].Vel,
                md.results.TransientSolution[11].Thickness,
                md.results.TransientSolution[11].SmbMassBalance,
                md.results.TransientSolution[11].FrictionWaterPressure,
                md.results.TransientSolution[23].Vel,
                md.results.TransientSolution[23].Thickness,
                md.results.TransientSolution[23].SmbMassBalance,
                md.results.TransientSolution[23].FrictionWaterPressure]
