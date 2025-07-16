import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class autodiff(class_registry.manage_state):
    '''
    autodiff Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.isautodiff = 0.
        self.dependents = 'List dependents'
        self.independents = 'List independents'
        self.driver = 'fos_forward'
        self.obufsize = np.nan
        self.lbufsize = np.nan
        self.cbufsize = np.nan
        self.tbufsize = np.nan
        self.gcTriggerMaxSize = np.nan
        self.gcTriggerRatio = np.nan
        self.tapeAlloc = np.nan
        self.outputTapeMemory = 0.
        self.outputTime = 0.
        self.enablePreaccumulation = 0.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):  # {{{
        s = '      automatic differentiation parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'isautodiff', "indicates if the automatic differentiation is activated"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'dependents', "list of dependent variables"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'independents', "list of independent variables"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'driver', "ADOLC driver ('fos_forward' or 'fov_forward')"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'obufsize', "Number of operations per buffer (== OBUFSIZE in usrparms.h)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'lbufsize', "Number of locations per buffer (== LBUFSIZE in usrparms.h)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'cbufsize', "Number of values per buffer (== CBUFSIZE in usrparms.h)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'tbufsize', "Number of taylors per buffer (<=TBUFSIZE in usrparms.h)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'gcTriggerRatio', "free location block sorting / consolidation triggered if the ratio between allocated and used locations exceeds gcTriggerRatio"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'gcTriggerMaxSize', "free location block sorting / consolidation triggered if the allocated locations exceed gcTriggerMaxSize)"))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'tapeAlloc', 'Iteration count of a priori memory allocation of the AD tape'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'outputTapeMemory', 'Write AD tape memory statistics to file ad_mem.dat'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'outputTime', 'Write AD recording and evaluation times to file ad_time.dat'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'enablePreaccumulation', 'Enable CoDiPack preaccumulation in augmented places'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - autodiff Class'
        return s
