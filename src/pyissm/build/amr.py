from . import build_utils
from . import class_registry

@class_registry.register_class
class amr(class_registry.manage_state):
    """
    Adaptive Mesh Refinement (AMR) class for ISSM.

    This class encapsulates parameters and configuration for adaptive mesh refinement (AMR) in the ISSM (Ice Sheet System Model) framework.
    It allows users to control mesh resolution based on various criteria such as velocity, thickness error, and deviatoric stress error, 
    as well as to refine mesh near grounding lines and ice fronts.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    hmin : float, default=100
        Minimum element length.
    hmax : float, default=100000.0
        Maximum element length.
    fieldname : str, default='Vel'
        Name of the input field used to compute the metric (should be an input of FemModel).
    err : int, default=3
        Error estimator type or flag.
    keepmetric : int, default=1
        Indicates whether the metric should be kept at every remeshing time.
    gradation : float, default=1.5
        Maximum ratio between two adjacent edge lengths.
    groundingline_resolution : float, default=500
        Element length near the grounding line.
    groundingline_distance : float, default=0
        Distance around the grounding line for mesh refinement.
    icefront_resolution : float, default=500
        Element length near the ice front.
    icefront_distance : float, default=0
        Distance around the ice front for mesh refinement.
    thicknesserror_resolution : float, default=500
        Element length when thickness error estimator is used.
    thicknesserror_threshold : float, default=0
        Maximum threshold for thickness error.
    thicknesserror_groupthreshold : float, default=0
        Maximum group threshold for thickness error.
    thicknesserror_maximum : float, default=0
        Maximum permitted thickness error.
    deviatoricerror_resolution : float, default=500
        Element length when deviatoric stress error estimator is used.
    deviatoricerror_threshold : float, default=0
        Maximum threshold for deviatoric stress error.
    deviatoricerror_groupthreshold : float, default=0
        Maximum group threshold for deviatoric stress error.
    deviatoricerror_maximum : float, default=0
        Maximum permitted deviatoric stress error.
    restart : int, default=0
        Indicates if ReMesh() should be called before the first time step.

    Methods
    -------
    __init__(self, other=None)
        Initializes the AMR parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the AMR parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    -------
    md.amr = pyissm.build.amr()
    md.amr.hmin = 50
    md.amr.fieldname = 'Thickness'
    """ 
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

