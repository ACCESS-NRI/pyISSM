#Test Name: SquareShelfDamageEvolutionSSA2dPralong
import numpy as np
import pyissm

# Parameterise model
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/Square.exp', 150000)
md = pyissm.model.param.set_mask(md, 'all', None)
md.materials = pyissm.model.classes.materials.damageice()
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareShelf.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')
md.cluster.np = 3
md.damage.isdamage = 1
md.damage.D = 0.1 * np.ones(md.mesh.numberofvertices)
md.damage.spcdamage = np.nan * np.ones(md.mesh.numberofvertices)
md.damage.law = 1

md.damage.c1 = 1.e-11
md.damage.c2 = 0.4
md.damage.c3 = 1.e-3
md.damage.healing = 0.4
md.damage.stress_threshold = 1.e5
md.damage.stabilization = 1

md.damage.requested_outputs = ['default', 'DamageF']

# Execute model
md = pyissm.model.execute.solve(md, 'DamageEvolution')

# Fields and tolerances to track changes
field_names = ['D', 'F']
field_tolerances = [1.e-13, 1.e-13]
field_values = [md.results.DamageEvolutionSolution.DamageDbar,
                md.results.DamageEvolutionSolution.DamageF]
