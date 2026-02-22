#######################################################
###### PINE ISLAND GLACIER PARAMETERIZATION FILE ######
#######################################################

import pyissm
import numpy as np

try:
    ## BY DEFAULT: Load data using the CCD
    import datapool as dp
    print(f"All required datasets will be loaded from the ACCESS Community Cryosphere Datapool on {pyissm.tools.config.get_hostname()}.")

    catalog = dp.catalog.DataCatalog()

    print(' -- LOADING DATASETS -- ')
    bm = catalog.load_dataset('measures_bedmachine_antarctica', version = 'v3')
    racmo = catalog.load_dataset('racmo2.3p2_monthly_27km_1979-2022')
    velocity = catalog.load_dataset('measures_insar_based_antarctica_ice_velocity_map', version = 'v2')
    aq1 = catalog.load_dataset('antarctic_geothermal_heat_flow_model_aq1', resolution = '20km')

    racmo_temp = racmo.copy()
    racmo_smb = racmo.copy()
    

except:
    ## WORKING LOCALLY: Load data directly using xarray
    print(f"The ACCESS Community Cryosphere Datapool is not accessible on {pyissm.tools.config.get_hostname()}. Please update the paths to the required input datasets in the pig_param.py parameterization file.")
    import xarray as xr

    ## To load datasets locally using xarray, please update the following paths
    velocity = xr.open_dataset("<PATH_TO_MEASURES_V2_ICE_VELOCITY_DATASET>")
    bm = xr.open_dataset("<PATH_TO_MEASURES_BEDMACHINCE_V3_DATASET>")
    racmo_temp = xr.open_dataset("<PATH_TO_RACMO2.3p3_MONTHLY_27KM_1979-2022_TSKIN_DATASET>")
    racmo_smb = xr.open_dataset("<PATH_TO_RACMO2.3p3_MONTHLY_27KM_1979-2022_SMB_DATASET>")
    aq1 = xr.open_dataset("<PATH_TO_AQ1_20KM_DATASET>")

## -------------- BEGIN PARAMETERIZATION HERE -------------- ##

# Parameters to change
friction_coefficient = 10
temp_change = 0

# Name and coordinate system
md.miscellaneous.name = 'PIG'
md.mesh.epsg = 3031

## ------ GEOMETRY ------ ##
print(' -- ASSIGNING GEOMETRY -- ')
md.geometry.base = pyissm.data.interp.xr_to_mesh(bm, 'bed', md.mesh.x, md.mesh.y)
md.geometry.surface = pyissm.data.interp.xr_to_mesh(bm, 'surface', md.mesh.x, md.mesh.y)
md.geometry.thickness = md.geometry.surface - md.geometry.base

# Use hydrostatic equilibrium on ice shelf
di = md.materials.rho_ice / md.materials.rho_water

## Get floating nodes
floating = md.mask.ocean_levelset < 0

## Apply floatation criterion on floating nodes and redefine base/thickness accordingly
md.geometry.thickness[floating] = 1 / (1 - di) * md.geometry.surface[floating]
md.geometry.base[floating] = md.geometry.surface[floating] - md.geometry.thickness[floating]

# Set minimum thickness to 1 m
md.geometry.thickness[md.geometry.thickness < 1] = 1
md.geometry.surface = md.geometry.thickness + md.geometry.base

# Set bed = base and lower bathymetry by 1000 m
md.geometry.bed = md.geometry.base.copy()
md.geometry.bed[floating]  = md.geometry.bed[floating] - 1000

## ------ INITIALIZATION ------ ##
print(' -- ASSIGNING INITIALIZATION FIELDS -- ')

# Interpolate racmo surface temperature
x, y = pyissm.tools.general.ll_to_xy(racmo_temp.lat, racmo_temp.lon, -1)
temp = racmo_temp['tskin']
temp = temp.sel(time = temp['time.year'] == 1995)
temp_ann = temp.mean('time').squeeze().values
temp_mesh = pyissm.data.interp.points_to_mesh(x, y, temp_ann, md.mesh.x, md.mesh.y)
md.initialization.temperature = np.minimum(temp_mesh, 273.15)  # Cap maximum temperature to 0 degC

# Interpolate surface velocities
md.initialization.vx = pyissm.data.interp.xr_to_mesh(velocity, 'VX', md.mesh.x, md.mesh.y, fill_nan = True)
md.initialization.vy = pyissm.data.interp.xr_to_mesh(velocity, 'VY', md.mesh.x, md.mesh.y, fill_nan = True)
md.initialization.vz = np.full(md.mesh.numberofvertices, 0)

# Add velocities to inversion class
md.inversion.vx_obs = md.initialization.vx.copy()
md.inversion.vy_obs = md.initialization.vy.copy()
md.inversion.vel_obs = np.sqrt(md.inversion.vx_obs ** 2 + md.inversion.vy_obs **2)

# Initialize pressure
md.initialization.pressure = md.materials.rho_ice * md.constants.g * md.geometry.thickness

# Contrust ice rheology properties
md.materials.rheology_n = 3 * np.ones(md.mesh.numberofelements, )
md.materials.rheology_B = pyissm.tools.materials.paterson(md.initialization.temperature)

## ------ FORCINGS ------ ##
print(' -- ASSIGNING FORCINGS -- ')

## SMB
x, y = pyissm.tools.general.ll_to_xy(racmo_smb.lat, racmo_smb.lon, -1)
smb = racmo_smb['smb']
smb = smb.sel(time = smb['time.year'] == 1995)
smb_ann = smb.groupby('time.year').sum('time').squeeze().values
smb_ann = smb_ann / md.materials.rho_ice  # Convert from m w.e. to m ice equivalent
# Interpolate to mesh vertices
smb_mesh = pyissm.data.interp.points_to_mesh(x, y, smb_ann, md.mesh.x, md.mesh.y)
md.smb.mass_balance = smb_mesh

## Geothermal heat flux
md.basalforcings.geothermalflux = pyissm.data.interp.xr_to_mesh(aq1,
                                                                'Q',
                                                                md.mesh.x,
                                                                md.mesh.y,
                                                                x_var = 'X',
                                                                y_var = 'Y',
                                                                crop_buffer = 20000, # Increase crop buffer due to coarse data grid
                                                                fill_nan = True) 

## Friction
md.friction.coefficient = friction_coefficient * np.ones(md.mesh.numberofvertices, )
md.friction.p = np.ones(md.mesh.numberofelements, )
md.friction.q = np.ones(md.mesh.numberofelements, )

## No friction on floating ice
md.friction.coefficient[md.mask.ocean_levelset < 0] = 0

## Enable sub-element grounding line migration
md.groundingline.migration = 'SubelementMigration'

## ------ BOUNDARY CONDITIONS ------ ##
print(' -- SETTING BOUNDARY CONDITIONS -- ')
md = pyissm.model.bc.set_marine_ice_sheet_bc(md)
md.basalforcings.floatingice_melting_rate = np.zeros(md.mesh.numberofvertices, )
md.basalforcings.groundedice_melting_rate = np.zeros(md.mesh.numberofvertices, )
md.thermal.spctemperature = md.initialization.temperature;
md.masstransport.spcthickness = np.full(md.mesh.numberofvertices, np.nan)

