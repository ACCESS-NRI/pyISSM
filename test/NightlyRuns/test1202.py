#Test Name: EISMINTStress1
import pyissm


"""
Test on the stressbalance model and the masstransport in 2d
"""

printingflag = False

#tests 3 and 4: using Glen's flow law
md = pyissm.model.Model()
md = pyissm.model.mesh.triangle(md, '../assets/Exp/SquareEISMINT.exp', 3550.)
md = pyissm.model.param.set_mask(md, 'all', None)
md = pyissm.model.param.parameterize(md, '../assets/Par/SquareEISMINT.py')
md = pyissm.model.param.set_flow_equation(md, SSA = 'all')  #SSA's model and 2d
#Compute solution for SSA's model
md.cluster.np = 8
md = pyissm.model.execute.solve(md, 'Stressbalance')

#plot results
vx = md.results.StressbalanceSolution.Vx
vy = md.results.StressbalanceSolution.Vy

#plotmodel(md, 'data', vx, 'contourlevels', {0, 20, 40, 60, 60, 100, 120, 140, 160, 180, -20, -40, -60, -80, -100, -120, -140, -160, -180}, ...
#       'contourcolor', 'k')
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('eismintdiag1vx', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 2, 'hardcopy', 'off')
#       system(['mv eismintdiag1vx.png ' ISSM_DIR '/website/doc_pdf/validation/Images/EISMINT/IceShelf '])

#plotmodel(md, 'data', vy, 'contourlevels', { -100, -200, -300, -400, -500, -600, -700, -800, -900, -1000}, ...
#       'contourcolor', 'k')
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('eismintdiag1vy', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 2, 'hardcopy', 'off')
#       system(['mv eismintdiag1vy.png ' ISSM_DIR '/website/doc_pdf/validation/Images/EISMINT/IceShelf '])

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy']
field_tolerances = [1e-13, 1e-13]
field_values = [vx, vy]
