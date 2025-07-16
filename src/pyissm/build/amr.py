from . import build_utils
from . import class_registry

@class_registry.register_class
class amr(class_registry.manage_state):
    '''
    amr Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.hmin = 100
        self.hmax = 100e3
        self.fieldname = 'Vel'
        self.err = 3
        self.keepmetric = 1
        self.gradation = 1.5
        self.groundingline_resolution = 500
        self.groundingline_distance = 0
        self.icefront_resolution = 500
        self.icefront_distance = 0
        self.thicknesserror_resolution = 500
        self.thicknesserror_threshold = 0
        self.thicknesserror_groupthreshold = 0
        self.thicknesserror_maximum = 0
        self.deviatoricerror_resolution = 500
        self.deviatoricerror_threshold = 0
        self.deviatoricerror_groupthreshold = 0
        self.deviatoricerror_maximum = 0
        self.restart = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   amr parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'hmin', 'minimum element length'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'hmax', 'maximum element length'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'fieldname', 'name of input that will be used to compute the metric (should be an input of FemModel)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'keepmetric', 'indicates whether the metric should be kept every remeshing time'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'gradation', 'maximum ratio between two adjacent edges'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'groundingline_resolution', 'element length near the grounding line'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'groundingline_distance', 'distance around the grounding line which elements will be refined'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'icefront_resolution', 'element length near the ice front'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'icefront_distance', 'distance around the ice front which elements will be refined'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thicknesserror_resolution', 'element length when thickness error estimator is used'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thicknesserror_threshold', 'maximum threshold thickness error permitted'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thicknesserror_groupthreshold', 'maximum group threshold thickness error permitted'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'thicknesserror_maximum', 'maximum thickness error permitted'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'deviatoricerror_resolution', 'element length when deviatoric stress error estimator is used'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'deviatoricerror_threshold', 'maximum threshold deviatoricstress error permitted'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'deviatoricerror_groupthreshold', 'maximum group threshold deviatoric stress error permitted'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'deviatoricerror_maximum', 'maximum deviatoricstress error permitted'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'restart', 'indicates if ReMesh() will call before first time step'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - amr Class'
        return s

