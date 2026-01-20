#Test Name: SquareShelfConstrainedStressSSA2d
import pyissm
import numpy as np

md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 50000)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelfConstrained.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3
#outputs
md.stressbalance.requested_outputs = ['default', 'DeviatoricStressxx', 'DeviatoricStressyy', 'DeviatoricStressxy', 'Outputdefinition1', 'Outputdefinition2', 'Outputdefinition3', 'Outputdefinition4', 'Outputdefinition5', 'Outputdefinition6']
md.outputdefinition.definitions = [pyissm.model.classes.massfluxatgate(name  = 'MassFlux1', profilename  = '../assets/Exp/MassFlux1.exp', definitionstring  ='Outputdefinition1'),
                                   pyissm.model.classes.massfluxatgate(name  = 'MassFlux2', profilename  = '../assets/Exp/MassFlux2.exp', definitionstring  ='Outputdefinition2'),
                                   pyissm.model.classes.massfluxatgate(name  = 'MassFlux3', profilename  = '../assets/Exp/MassFlux3.exp', definitionstring  ='Outputdefinition3'),
                                   pyissm.model.classes.massfluxatgate(name  = 'MassFlux4', profilename  = '../assets/Exp/MassFlux4.exp', definitionstring  ='Outputdefinition4'),
                                   pyissm.model.classes.massfluxatgate(name  = 'MassFlux5', profilename  = '../assets/Exp/MassFlux5.exp', definitionstring  ='Outputdefinition5'),
                                   pyissm.model.classes.massfluxatgate(name  = 'MassFlux6', profilename  = '../assets/Exp/MassFlux6.exp', definitionstring  ='Outputdefinition6')]
for d in md.outputdefinition.definitions:
    print(d.name, "segments is nan?", np.isnan(d.segments).all() if hasattr(d, "segments") else None)

md = pyissm.model.execute.solve(md, 'Stressbalance')
sol = md.results.StressbalanceSolution
print([a for a in dir(sol) if "Outputdefinition" in a or "MassFlux" in a])
#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vel', 'Pressure',
               'DeviatoricStressxx', 'DeviatoricStressyy', 'DeviatoricStressxy',
               'MassFlux1', 'MassFlux2', 'MassFlux3', 'MassFlux4', 'MassFlux5', 'MassFlux6']
field_tolerances = [4e-13, 4e-13, 4e-13, 1e-13,
                    2e-13, 2e-13, 2e-13,
                    1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure,
                md.results.StressbalanceSolution.DeviatoricStressxx,
                md.results.StressbalanceSolution.DeviatoricStressyy,
                md.results.StressbalanceSolution.DeviatoricStressxy,
                md.results.StressbalanceSolution.MassFlux1,
                md.results.StressbalanceSolution.MassFlux2,
                md.results.StressbalanceSolution.MassFlux3,
                md.results.StressbalanceSolution.MassFlux4,
                md.results.StressbalanceSolution.MassFlux5,
                md.results.StressbalanceSolution.MassFlux6]
