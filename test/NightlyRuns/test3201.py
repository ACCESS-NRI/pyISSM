# Test Name: SquareShelfTransientCalibrationNBEcodipack
import numpy as np
import pyissm

# -----------------------------
# Generate observations
# -----------------------------
md = pyissm.model.Model()
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../Exp/Square.exp', 100000.)
md = pyissm.model.param.set_mask(md, 'all', None)  # setmask(md,'all','')
md = pyissm.model.param.parameterize(md, '../Par/SquareShelf.par')
md = pyissm.model.param.set_flow_equation(md, SSA='all')
md.cluster = pyissm.cluster.generic('np', 2)

# -----------------------------
# Create "true" time series for B
# -----------------------------
md.timestepping.interp_forcing = 0
md.timestepping.final_time = 2.0 * md.timestepping.time_step

# rheology_B is (nelems, 2) time series values, then an extra "time row" [0.01, 2*dt]
ne = md.mesh.numberofelements
B = 1.8e8 * np.ones((ne, 2))

# MATLAB: mean(md.mesh.x(md.mesh.elements),2) < mean(md.mesh.y(md.mesh.elements),2)  (on elements)
# Note: ISSM elements are 1-based in MATLAB; in pyISSM they are usually 1-based too.
# We'll convert to 0-based indices for numpy indexing.
ele = md.mesh.elements.copy()
if ele.min() == 1:
    ele0 = ele - 1
else:
    ele0 = ele

x_e_mean = md.mesh.x[ele0].mean(axis=1)
y_e_mean = md.mesh.y[ele0].mean(axis=1)
mask = np.where(x_e_mean < y_e_mean)[0]
B[mask, 1] = 1.4e8

# append time row
dt = md.timestepping.time_step
B = np.vstack([B, np.array([0.01, 2.0 * dt])])

md.materials.rheology_B = B

# -----------------------------
# Initial values
# -----------------------------
nv = md.mesh.numberofvertices
md.initialization.vx = np.zeros((nv,))
md.initialization.vy = np.zeros((nv,))
md.initialization.pressure = np.zeros((nv,))
md.initialization.temperature = np.zeros((nv,))

md.basalforcings.geothermalflux = np.zeros((nv,))
md.thermal.spctemperature = np.full((nv,), np.nan)

# solve forward transient to create observations
md = pyissm.execute.solve(md, 'tr')

# -----------------------------
# Modify rheology, now constant
# -----------------------------
md.materials.rheology_B[:-1, :] = 1.8e8

# -----------------------------
# Set cost functions (one set per transient time)
# -----------------------------
md.outputdefinition.definitions = []
md.autodiff.dependents = []

count = 1
for sol in md.results.TransientSolution:
    vx_obs = sol.Vx
    vy_obs = sol.Vy
    z_obs  = sol.Surface
    time   = sol.time
    weights = np.ones((nv,))

    # 1) LogVel misfit
    od1 = pyissm.model.classes.cfsurfacelogvel()
    od1.name = f'LogVelMis{count}'
    od1.definitionstring = f'Outputdefinition{count}'
    od1.vxobs_string = 'VxObs'
    od1.vxobs = vx_obs
    od1.vyobs_string = 'VyObs'
    od1.vyobs = vy_obs
    od1.weights = weights
    od1.wei
