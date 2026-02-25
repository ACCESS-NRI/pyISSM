#Test Name: SquareSheetShelfTranMeltFCT 
import pyissm
import copy
import numpy as np

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000.)
md = pyissm.model.param.set_mask(md, '../assets/Exp/SquareShelf.exp', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareSheetShelf.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3
md.geometry.bed = copy.deepcopy(md.geometry.base)
pos = np.nonzero(md.mask.ocean_levelset < 0.)
md.geometry.bed[pos] = md.geometry.bed[pos] - 10
md.friction = pyissm.model.classes.friction.coulomb()
md.friction.coefficient = 20 * np.ones(md.mesh.numberofvertices)
md.friction.p = 1 * np.ones(md.mesh.numberofelements)
md.friction.q = 1 * np.ones(md.mesh.numberofelements)
md.friction.coefficientcoulomb = 0.02 * np.ones(md.mesh.numberofvertices)
md.transient.isthermal = False
md.transient.isgroundingline = True
md.transient.requested_outputs = ['default', 'GroundedArea', 'FloatingArea', 'TotalFloatingBmb', 'TotalGroundedBmb', 'TotalSmb']
md.masstransport.stabilization = 4

# Execute model
md = pyissm.model.execute.solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = [
    'Vx1', 'Vy1', 'Vel1', 'Pressure1', 'Bed1', 'Surface1', 'Thickness1', 'GroundedArea1', 'FloatingArea1', 'TotalFloatingBmb1', 'TotalGroundedBmb1', 'TotalSmb1', 
    'Vx2', 'Vy2', 'Vel2', 'Pressure2', 'Bed2', 'Surface2', 'Thickness2', 'GroundedArea2', 'FloatingArea2', 'TotalFloatingBmb2', 'TotalGroundedBmb2', 'TotalSmb2', 
    'Vx3', 'Vy3', 'Vel3', 'Pressure3', 'Bed3', 'Surface3', 'Thickness3', 'GroundedArea3', 'FloatingArea3', 'TotalFloatingBmb3', 'TotalGroundedBmb3', 'TotalSmb3'
]
field_tolerances = [
    2e-13, 2e-13, 2e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 
    2e-13, 2e-13, 2e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 
    2e-13, 2e-13, 2e-13, 1e-13, 2e-13, 5e-13, 1e-12, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 
]
field_values = [
    md.results.TransientSolution[0].Vx,
    md.results.TransientSolution[0].Vy,
    md.results.TransientSolution[0].Vel,
    md.results.TransientSolution[0].Pressure,
    md.results.TransientSolution[0].Base,
    md.results.TransientSolution[0].Surface,
    md.results.TransientSolution[0].Thickness,
    md.results.TransientSolution[0].GroundedArea,
    md.results.TransientSolution[0].FloatingArea,
    md.results.TransientSolution[0].TotalFloatingBmb,
    md.results.TransientSolution[0].TotalGroundedBmb,
    md.results.TransientSolution[0].TotalSmb,
    md.results.TransientSolution[1].Vx,
    md.results.TransientSolution[1].Vy,
    md.results.TransientSolution[1].Vel,
    md.results.TransientSolution[1].Pressure,
    md.results.TransientSolution[1].Base,
    md.results.TransientSolution[1].Surface,
    md.results.TransientSolution[1].Thickness,
    md.results.TransientSolution[1].GroundedArea,
    md.results.TransientSolution[1].FloatingArea,
    md.results.TransientSolution[1].TotalFloatingBmb,
    md.results.TransientSolution[1].TotalGroundedBmb,
    md.results.TransientSolution[1].TotalSmb,
    md.results.TransientSolution[2].Vx,
    md.results.TransientSolution[2].Vy,
    md.results.TransientSolution[2].Vel,
    md.results.TransientSolution[2].Pressure,
    md.results.TransientSolution[2].Base,
    md.results.TransientSolution[2].Surface,
    md.results.TransientSolution[2].Thickness,
    md.results.TransientSolution[2].GroundedArea,
    md.results.TransientSolution[2].FloatingArea,
    md.results.TransientSolution[2].TotalFloatingBmb,
    md.results.TransientSolution[2].TotalGroundedBmb,
    md.results.TransientSolution[2].TotalSmb
]
