#Test Name: GiaIvinsBenchmarksAB2dC
from socket import gethostname

import numpy as np

import pyissm

# Benchmark experiments (Figure A2a Ivins and James, 1999, Geophys. J. Int.)
md = pyissm.model.mesh.triangle(pyissm.model.Model(), '../assets/Exp/RoundFrontEISMINT.exp', 200000.)
md = pyissm.model.param.set_mask(md, None, None)
md = pyissm.model.param.parameterize(md, '../assets/Par/GiaIvinsBenchmarksCD.py')

# indicate what you want to compute
md.gia.cross_section_shape = 1  # for square-edged x-section

# evaluation time (termed start_time)
md.timestepping.start_time = 0.3  # for t \approx 0 kyr : to get elastic response!
md.timestepping.final_time = 2500000  # 2,500 kyr

# define loading history
md.geometry.thickness = np.array([
    np.append(md.geometry.thickness * 0.0, 0.0),
    np.append(md.geometry.thickness / 2.0, 0.1),
    np.append(md.geometry.thickness, 0.2),
    np.append(md.geometry.thickness, md.timestepping.start_time)
    ]).T

# find out elements that have zero loads throughout the loading history
pos = np.where(np.abs(md.geometry.thickness[0:-2, :].sum(axis=1)) == 0)[0]
md.mask.ice_levelset[pos] = 1 # no ice

md.cluster.np = 3

# solve for GIA deflection
md = pyissm.model.execute.solve(md, 'Gia')

# Test Name: GiaIvinsBenchmarksAB2dC1
U_AB2dC1 = md.results.GiaSolution.UGia
URate_AB2dC1 = md.results.GiaSolution.UGiaRate

# Test Name: GiaIvinsBenchmarksAB2dC2
# different evaluation time # {{{
md.timestepping.start_time = 1000.3 # for t \approx 1 kyr
md.geometry.thickness[-1, -1] = md.timestepping.start_time

md = pyissm.model.execute.solve(md, 'Gia')

U_AB2dC2 = md.results.GiaSolution.UGia
URate_AB2dC2 = md.results.GiaSolution.UGiaRate
# }}}

# Test Name: GiaIvinsBenchmarksAB2dC3
# different evaluation time # {{{
md.timestepping.start_time = 2400000 # for t \approx \infty
md.geometry.thickness[-1, -1] = md.timestepping.start_time

md = pyissm.model.execute.solve(md, 'Gia')

U_AB2dC3 = md.results.GiaSolution.UGia
URate_AB2dC3 = md.results.GiaSolution.UGiaRate
# }}}

# Fields and tolerances to track changes
field_names = ['U_AB2dC1', 'URate_AB2dC1', 'U_AB2dC2', 'URate_AB2dC2', 'U_AB2dC3', 'URate_AB2dC3']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [U_AB2dC1, URate_AB2dC1, U_AB2dC2, URate_AB2dC2, U_AB2dC3, URate_AB2dC3]
