#Test Name: SquareShelfConstrainedStressMOLHO2d
import pyissm

md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 50000)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
md = pyissm.model.param.set_flow_equation(md, MOLHO = 'all')
md.cluster.np = 2
#outputs
#FIXME compute the stress components for MOLHO
md.stressbalance.requested_outputs = ['default', 'VxSurface', 'VySurface', 'VxShear', 'VyShear', 'VxBase', 'VyBase', 'MassFlux1', 'MassFlux2', 'MassFlux3', 'MassFlux4', 'MassFlux5', 'MassFlux6']
#md.stressbalance.requested_outputs = ['default', 'DeviatoricStressxx', 'DeviatoricStressyy', 'DeviatoricStressxy', 'MassFlux1', 'MassFlux2', 'MassFlux3', 'MassFlux4', 'MassFlux5', 'MassFlux6']
md.outputdefinition.definitions = [pyissm.model.classes.massfluxatgate(name="MassFlux1", profilename="../assets/Exp/MassFlux1.exp", definitionstring="Outputdefinition1"),
                                   pyissm.model.classes.massfluxatgate(name="MassFlux2", profilename="../assets/Exp/MassFlux2.exp", definitionstring="Outputdefinition2"),
                                   pyissm.model.classes.massfluxatgate(name="MassFlux3", profilename="../assets/Exp/MassFlux3.exp", definitionstring="Outputdefinition3"),
                                   pyissm.model.classes.massfluxatgate(name="MassFlux4", profilename="../assets/Exp/MassFlux4.exp", definitionstring="Outputdefinition4"),
                                   pyissm.model.classes.massfluxatgate(name="MassFlux5", profilename="../assets/Exp/MassFlux5.exp", definitionstring="Outputdefinition5"),
                                   pyissm.model.classes.massfluxatgate(name="MassFlux6", profilename="../assets/Exp/MassFlux6.exp", definitionstring="Outputdefinition6")]
md = pyissm.model.bc.set_molho_bc(md)
md = pyissm.model.execute.solve(md, 'Stressbalance')

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vel', 'Pressure', 'VxSurface', 'VySurface', 'VxShear', 'VyShear', 'VxBase', 'VyBase',
               'MassFlux1', 'MassFlux2', 'MassFlux3', 'MassFlux4', 'MassFlux5', 'MassFlux6']
#'DeviatoricStressxx', 'DeviatoricStressyy', 'DeviatoricStressxy',
field_tolerances = [5e-13, 6e-13, 6e-13, 1e-13, 5e-13, 6e-13, 2e-13, 4e-13, 5e-13, 6e-13,
                    2e-13, 1e-13, 2e-13,
                    1e-13, 1e-13, 1e-13]
                    #1e-13, 1e-13, 1e-13]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure,
                md.results.StressbalanceSolution.VxSurface,
                md.results.StressbalanceSolution.VySurface,
                md.results.StressbalanceSolution.VxShear,
                md.results.StressbalanceSolution.VyShear,
                md.results.StressbalanceSolution.VxBase,
                md.results.StressbalanceSolution.VyBase,
                md.results.StressbalanceSolution.MassFlux1,
                md.results.StressbalanceSolution.MassFlux2,
                md.results.StressbalanceSolution.MassFlux3,
                md.results.StressbalanceSolution.MassFlux4,
                md.results.StressbalanceSolution.MassFlux5,
                md.results.StressbalanceSolution.MassFlux6]
#md.results.StressbalanceSolution.DeviatoricStressxx,
#md.results.StressbalanceSolution.DeviatoricStressyy,
#md.results.StressbalanceSolution.DeviatoricStressxy,
