#Test Name: PigTranARMAandStochasticforcings
import numpy as np
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Pig.exp', 8000)
md = pyissm.model.param.set_mask(md, '../assets/Exp/PigShelves.exp', '../assets/Exp/PigIslands.exp')
md = pyissm.model.param.parameterize(md, '../assets/Par/Pig.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.timestepping.start_time = 0
md.timestepping.time_step = 1
md.timestepping.final_time = 10
md.cluster.np = 2

# Basin separation
idb = np.zeros((md.mesh.numberofelements,))
iid1 = np.where(md.mesh.x >= -1.6e6)[0]
for ii in range(md.mesh.numberofelements):
    for vertex in range(3):
        if md.mesh.elements[ii][vertex] - 1 in iid1:  # one vertex in basin 1; NOTE: offset because of 1-based vertex indexing
            idb[ii] = 1
    if idb[ii] == 0:  # no vertex was found in basin 1
        for vertex in range(3):
            idb[ii] = 2
nb_bas = 2

# SMB
numparams  = 1
numbreaks  = 0
intercept     = np.array([[0.5],[0.01]])
polynomialparams = np.copy(intercept)
datebreaks = np.nan
md.smb = pyissm.model.classes.smb.arma()
md.smb.num_basins = nb_bas  # number of basins
md.smb.basin_id = idb  # prescribe basin ID number to elements;
md.smb.num_params = 1 * numparams
md.smb.num_breaks = 1 * numbreaks
md.smb.polynomialparams = 1 * polynomialparams
md.smb.datebreaks = 1 * datebreaks
md.smb.ar_order = 4
md.smb.ma_order = 4
md.smb.arma_timestep = 2.0  #timestep of the ARMA model [yr]
md.smb.arlag_coefs = np.array([[0.2,0.1,0.05,0.01],[0.4,0.2,-0.2,0.1]])
md.smb.malag_coefs = np.array([[0.1,0.1,0.2,0.3],[0.5,0.8,1.3,2.4]])

# Calving
md.mask.ice_levelset = 1e4*(md.mask.ice_levelset + 0.5)
md.calving.calvingrate = 0.1 * np.ones((md.mesh.numberofvertices,))
md.levelset.spclevelset = np.full((md.mesh.numberofvertices,), np.nan)
md.levelset.migration_max = 10.0
md.frontalforcings.meltingrate = np.zeros((md.mesh.numberofvertices,))

# Basal forcing implementation
numparams = 2
numbreaks = 1
intercept = np.array([[3.0,4.0],[1.0,0.5]])
trendlin  = np.array([[0.0,0.1],[0,0]])
polynomialparams = np.transpose(np.stack((intercept,trendlin)),(1,2,0))
datebreaks = np.array([[6],[7]])

md.basalforcings = pyissm.model.classes.basalforcings.lineararma()
md.basalforcings.num_basins = nb_bas
md.basalforcings.basin_id  = idb
md.basalforcings.const = np.array([[1.0, 2.50]])  # intercept values of DeepwaterMelt in basins [m/yr]
md.basalforcings.trend  = np.array([[0.2, 0.01]])  # trend values of DeepwaterMelt in basins [m/yr^2]
md.basalforcings.arma_initialtime = md.timestepping.start_time  # initial time in the ARMA model parameterization [yr]
md.basalforcings.ar_order = 1
md.basalforcings.ma_order = 1
md.basalforcings.polynomialparams = 1 * polynomialparams
md.basalforcings.datebreaks = 1 * datebreaks
md.basalforcings.num_params = 1 * numparams
md.basalforcings.num_breaks = 1 * numbreaks
md.basalforcings.arma_timestep = 1.0  # timestep of the ARMA model [yr]
md.basalforcings.arlag_coefs = np.array([[0.0], [0.1]])  # autoregressive parameters
md.basalforcings.malag_coefs = np.array([[0.55], [0.34]])  # moving-average parameters
md.basalforcings.deepwater_elevation = np.array([[-1000, -1520]])
md.basalforcings.upperwater_elevation = np.array([[0, -50]])
md.basalforcings.upperwater_melting_rate = np.array([[0,0]])
md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices,))

# Covariance matrices
sdvsmb  = np.array([1,1])
sdvclv  = np.array([0.1,0.01])
sdvdwm  = np.array([300,300])
vecsdv  = np.concatenate((sdvsmb,sdvclv,sdvdwm))
corrmat = np.array([[1.0, 0., 0., 0., 0., 0.],
                    [0., 1.0, 0., 0., 0., 0.],
                    [0., 0., 1.0, 0.4, 0.1, 0.1],
                    [0., 0., 0.4, 1.0, 0.1, 0.1],
                    [0., 0., 0.1, 0.1, 1.0, 0.3],
                    [0., 0., 0.1, 0.1, 0.3, 1.0]])
covglob0 = np.diag(vecsdv) @ corrmat @ np.diag(vecsdv)
covglob1 = 2*covglob0
multcov  = np.stack((covglob0,covglob1),axis=2) 
tmcov    = np.array([[0,5]])

# Stochastic forcing
md.stochasticforcing.isstochasticforcing = 1
md.stochasticforcing.fields = ['SMBarma', 'DefaultCalving', 'BasalforcingsDeepwaterMeltingRatearma']
md.stochasticforcing.defaultdimension = 2
md.stochasticforcing.default_id = idb
md.stochasticforcing.covariance = multcov  # global covariance among- and between-fields
md.stochasticforcing.timecovariance = tmcov  #simulation times when covariance matrix switches 
md.stochasticforcing.randomflag = 0  # determines true/false randomness

md.transient.ismovingfront = 1
md.transient.requested_outputs = ['default', 'SmbMassBalance', 'BasalforcingsFloatingiceMeltingRate', 'BasalforcingsSpatialDeepwaterMeltingRate']
md.transient.isstressbalance = 1
md.transient.ismasstransport = 1
md.transient.issmb = 1
md.transient.isthermal = 0
md.transient.isgroundingline = 1

# Execute model
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = [
    'Vx1', 'Vy1', 'Vel1', 'Thickness1', 'SmbMassBalance1', 'BasalforcingsFloatingiceMeltingRate1', 'BasalforcingsSpatialDeepwaterMeltingRate1',
    'Vx5', 'Vy5', 'Vel5', 'Thickness5', 'SmbMassBalance5', 'BasalforcingsFloatingiceMeltingRate5', 'BasalforcingsSpatialDeepwaterMeltingRate5',
    'Vx10', 'Vy10', 'Vel10', 'Thickness10', 'SmbMassBalance10', 'BasalforcingsFloatingiceMeltingRate10', 'BasalforcingsSpatialDeepwaterMeltingRate10']

field_tolerances = [
    1e-11, 1e-11, 2e-11, 1e-11, 1e10, 1e-9, 1e-10,
    1e-11, 1e-11, 2e-11, 9e-11, 1e10, 1e-9, 1e-10,
    2e-10, 2e-10, 2e-10, 1e-10, 1e10, 1e-9, 1e-10]
field_values = [
    md.results.TransientSolution[0].Vx,
    md.results.TransientSolution[0].Vy,
    md.results.TransientSolution[0].Vel,
    md.results.TransientSolution[0].Thickness,
    md.results.TransientSolution[0].SmbMassBalance,
    md.results.TransientSolution[0].BasalforcingsFloatingiceMeltingRate,
    md.results.TransientSolution[0].BasalforcingsSpatialDeepwaterMeltingRate,
    md.results.TransientSolution[4].Vx,
    md.results.TransientSolution[4].Vy,
    md.results.TransientSolution[4].Vel,
    md.results.TransientSolution[4].Thickness,
    md.results.TransientSolution[4].SmbMassBalance,
    md.results.TransientSolution[4].BasalforcingsFloatingiceMeltingRate,
    md.results.TransientSolution[4].BasalforcingsSpatialDeepwaterMeltingRate,
    md.results.TransientSolution[9].Vx,
    md.results.TransientSolution[9].Vy,
    md.results.TransientSolution[9].Vel,
    md.results.TransientSolution[9].Thickness,
    md.results.TransientSolution[9].SmbMassBalance,
    md.results.TransientSolution[9].BasalforcingsFloatingiceMeltingRate,
    md.results.TransientSolution[9].BasalforcingsSpatialDeepwaterMeltingRate]
