#Test Name: SquareNoDynHydrologyDCSmbCoupled
import pyissm
import numpy as np

md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 100000.)
md = pyissm.model.param.set_mask(md, None, None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareNoDyn.py')
md.cluster.np = 1

md.transient.ishydrology = True
md.transient.issmb = True
md.hydrology = pyissm.model.classes.hydrology.dc()
md.hydrology = md.hydrology.initialise(md)
md.smb = pyissm.model.classes.smb.gradientscomponents()

md.hydrology.isefficientlayer = 1

md.hydrology.sedimentlimit_flag = 1
md.hydrology.sedimentlimit = 400.0
md.hydrology.sediment_transmitivity = 3.0 * np.ones((md.mesh.numberofvertices))
md.hydrology.mask_thawed_node = np.ones((md.mesh.numberofvertices))

md.hydrology.mask_eplactive_node = np.zeros((md.mesh.numberofvertices))
md.hydrology.epl_conductivity = 3.
md.hydrology.epl_initial_thickness = 20
md.hydrology.epl_colapse_thickness = 1.0e-3
md.hydrology.epl_thick_comp = 0
md.hydrology.epl_max_thickness = 1

md.hydrology.spcsediment_head = np.nan * np.ones((md.mesh.numberofvertices))
md.hydrology.spcepl_head = np.nan * np.ones((md.mesh.numberofvertices))

md.initialization.sediment_head = np.zeros((md.mesh.numberofvertices))
md.initialization.epl_head = np.zeros((md.mesh.numberofvertices))
md.initialization.epl_thickness = np.ones((md.mesh.numberofvertices))

#try to keep the different steps as common multipliers of eachother
#for the sake of precision the values of the model used as input should be on a shorter time-step-
#you can plot the results of this test and check how the time stepping is doing (see bellow)
md.hydrology.steps_per_step = 5
md.smb.steps_per_step = 10  #md.hydrology.steps_per_step
md.timestepping.time_step = 1.
md.timestepping.final_time = 20.0

smb_step = md.timestepping.time_step / md.smb.steps_per_step
hydro_step = md.timestepping.time_step / md.hydrology.steps_per_step
duration = np.arange(md.timestepping.start_time, md.timestepping.final_time + smb_step, smb_step)

ddf = 10.0e-3
md.smb.accuref = np.array([[0.5, 0.5], [md.timestepping.start_time, md.timestepping.final_time]])
md.smb.accualti = 0.0
md.smb.accugrad = np.array([[0., 0.], [md.timestepping.start_time, md.timestepping.final_time]])

#md.smb.runoffref = 9. * ddf * np.ones(np.shape(duration))  #constant input for testing
md.smb.runoffref = 0.9 * duration * ddf
md.smb.runoffref = np.vstack((md.smb.runoffref, duration))
md.smb.runoffalti = 0.0
md.smb.runoffgrad = np.array([[-6.5e-3 * ddf, -6.5e-3 * ddf], [md.timestepping.start_time, md.timestepping.final_time]])  #lapse rate *ddf*day per year

md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices))
md.basalforcings.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices))

# Execute model
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['SedimentWaterHead1', 'EplWaterHead1', 'SedimentHeadResidual1',
               'SedimentWaterHead4', 'EplWaterHead4', 'SedimentHeadResidual4',
               'SedimentWaterHead5', 'EplWaterHead5', 'SedimentHeadResidual5',
               'SedimentWaterHead9', 'EplWaterHead9', 'SedimentHeadResidual9',
               'EplWaterHead20', 'EplWaterHeadSubstep20', 'SedimentWaterHead20',
               'SedimentWaterHeadSubstep20']
field_tolerances = [1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13,
                    1e-13, 5e-12, 1e-11,
                    1e-13, 5e-12, 1e-11,
                    1e-13, 1e-13, 1e-13,
                    1e-13]
field_values = [md.results.TransientSolution[0].SedimentHead,
                md.results.TransientSolution[0].EplHead,
                md.results.TransientSolution[0].SedimentHeadResidual,
                md.results.TransientSolution[3].SedimentHead,
                md.results.TransientSolution[3].EplHead,
                md.results.TransientSolution[3].SedimentHeadResidual,
                md.results.TransientSolution[4].SedimentHead,
                md.results.TransientSolution[4].EplHead,
                md.results.TransientSolution[4].SedimentHeadResidual,
                md.results.TransientSolution[8].SedimentHead,
                md.results.TransientSolution[8].EplHead,
                md.results.TransientSolution[8].SedimentHeadResidual,
                md.results.TransientSolution[-1].EplHead,
                md.results.TransientSolution[-1].EplHeadSubstep,
                md.results.TransientSolution[-1].SedimentHead,
                md.results.TransientSolution[-1].SedimentHeadSubstep]