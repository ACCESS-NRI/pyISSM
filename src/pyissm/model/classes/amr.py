from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm.model import execute

@class_registry.register_class
class amr(class_registry.manage_state):
    """
    Adaptive Mesh Refinement (AMR) class for ISSM

    This class contains parameters and configuration options for adaptive mesh refinement (AMR) in the ISSM framework.
    It allows users to control mesh resolution based on various criteria such as velocity, thickness error,
    and deviatoric stress error, as well as to refine mesh near grounding lines and ice fronts.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default
        values, they will override the default values.

    Attributes
    ----------
    hmin : :class:`float`, default=100
        Minimum element length.
    hmax : :class:`float`, default=100000
        Maximum element length.
    fieldname : :class:`str`, default='Vel'
        Name of the input field used to compute the metric (should be an input of FemModel).
    err : :class:`int`, default=3
        Error estimator type or flag.
    keepmetric : :class:`int`, default=1
        Indicates whether the metric should be kept at every remeshing time.
    gradation : :class:`float`, default=1.5
        Maximum ratio between two adjacent edge lengths.
    groundingline_resolution : :class:`float`, default=500
        Element length near the grounding line.
    groundingline_distance : :class:`float`, default=0
        Distance around the grounding line for mesh refinement.
    icefront_resolution : :class:`float`, default=500
        Element length near the ice front.
    icefront_distance : :class:`float`, default=0
        Distance around the ice front for mesh refinement.
    thicknesserror_resolution : :class:`float`, default=500
        Element length when thickness error estimator is used.
    thicknesserror_threshold : :class:`float`, default=0
        Maximum threshold for thickness error.
    thicknesserror_groupthreshold : :class:`float`, default=0
        Maximum group threshold for thickness error.
    thicknesserror_maximum : :class:`float`, default=0
        Maximum permitted thickness error.
    deviatoricerror_resolution : :class:`float`, default=500
        Element length when deviatoric stress error estimator is used.
    deviatoricerror_threshold : :class:`float`, default=0
        Maximum threshold for deviatoric stress error.
    deviatoricerror_groupthreshold : :class:`float`, default=0
        Maximum group threshold for deviatoric stress error.
    deviatoricerror_maximum : :class:`float`, default=0
        Maximum permitted deviatoric stress error.
    restart : :class:`int`, default=0
        Indicates if ReMesh() should be called before the first time step.

    Methods
    -------
    __init__(self, other=None)
        Initializes the AMR parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the AMR parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    -------
    .. code-block:: python

        >>> md.amr = pyissm.model.classes.amr()
        >>> md.amr.hmin = 50
        >>> md.amr.fieldname = 'Thickness'
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

        s += '{}\n'.format(class_utils.fielddisplay(self, 'hmin', 'minimum element length'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'hmax', 'maximum element length'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'fieldname', 'name of input that will be used to compute the metric (should be an input of FemModel)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'keepmetric', 'indicates whether the metric should be kept every remeshing time'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'gradation', 'maximum ratio between two adjacent edges'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'groundingline_resolution', 'element length near the grounding line'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'groundingline_distance', 'distance around the grounding line which elements will be refined'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'icefront_resolution', 'element length near the ice front'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'icefront_distance', 'distance around the ice front which elements will be refined'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'thicknesserror_resolution', 'element length when thickness error estimator is used'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'thicknesserror_threshold', 'maximum threshold thickness error permitted'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'thicknesserror_groupthreshold', 'maximum group threshold thickness error permitted'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'thicknesserror_maximum', 'maximum thickness error permitted'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'deviatoricerror_resolution', 'element length when deviatoric stress error estimator is used'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'deviatoricerror_threshold', 'maximum threshold deviatoricstress error permitted'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'deviatoricerror_groupthreshold', 'maximum group threshold deviatoric stress error permitted'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'deviatoricerror_maximum', 'maximum deviatoricstress error permitted'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'restart', 'indicates if ReMesh() will call before first time step'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - amr Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [amr] parameters.
        
        Parameters
        ----------
        md : :class:`pyissm.model.Model`
            The model object to check.
        solution : :class:`pyissm.model.solution`
            The solution object to check.
        analyses : list of :class:`str`
            List of analyses to check consistency for.
        
        Returns 
        -------
        md : :class:`pyissm.model.Model`
            The model object with any consistency errors noted.
        """
        
        class_utils.check_field(md, fieldname = "amr.hmax", scalar = True, gt = 0, allow_nan = False)
        class_utils.check_field(md, fieldname = "amr.hmin", scalar = True, gt = 0, lt = self.hmax, allow_nan = False)
        class_utils.check_field(md, fieldname = "amr.keepmetric", scalar = True, ge = 0, le = 1, allow_nan = False)
        class_utils.check_field(md, fieldname = "amr.gradation", scalar = True, ge = 1.1, le = 5, allow_nan = False)
        class_utils.check_field(md, fieldname = "amr.groundingline_resolution", scalar = True, gt = 0, lt = self.hmax, allow_nan = False)
        class_utils.check_field(md, fieldname = "amr.groundingline_distance", scalar = True, ge = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "amr.icefront_resolution", scalar = True, gt = 0, lt = self.hmax, allow_nan = False)
        class_utils.check_field(md, fieldname = "amr.icefront_distance", scalar = True, ge = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "amr.thicknesserror_resolution", scalar = True, gt = 0, lt = self.hmax, allow_nan = False)
        class_utils.check_field(md, fieldname = "amr.thicknesserror_threshold", scalar = True, ge = 0, le = 1, allow_nan = False)
        class_utils.check_field(md, fieldname = "amr.thicknesserror_groupthreshold", scalar = True, ge = 0, le = 1, allow_nan = False)
        class_utils.check_field(md, fieldname = "amr.thicknesserror_maximum", scalar = True, ge = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "amr.deviatoricerror_resolution", scalar = True, gt = 0, lt = self.hmax, allow_nan = False)
        class_utils.check_field(md, fieldname = "amr.deviatoricerror_threshold", scalar = True, ge = 0, le = 1, allow_nan = False)
        class_utils.check_field(md, fieldname = "amr.deviatoricerror_groupthreshold", scalar = True, ge = 0, le = 1, allow_nan = False)
        class_utils.check_field(md, fieldname = "amr.deviatoricerror_maximum", scalar = True, ge = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "amr.restart", scalar = True, ge = 0, le = 1, allow_nan = False)

        return md

    # Marshall method for saving the amr parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [amr] parameters to a binary file.

        Parameters
        ----------
        fid : :class:`file object`
            The file object to write the binary data to.
        prefix : :class:`str`
            Prefix string used for data identification in the binary file.
        md : :class:`pyissm.model.Model`, optional
            ISSM model object needed in some cases.
            
        Returns
        -------
        None
        """

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.amr.type', data = 1, format = 'Integer')
        
        ## Write double fields
        fieldnames = ['hmin', 'hmax', 'err', 'gradation', 'groundingline_resolution', 'groundingline_distance', 'icefront_resolution',
                      'icefront_distance', 'thicknesserror_resolution', 'thicknesserror_threshold', 'thicknesserror_groupthreshold',
                      'thicknesserror_maximum', 'deviatoricerror_resolution', 'deviatoricerror_threshold', 'deviatoricerror_groupthreshold',
                      'deviatoricerror_maximum']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Double')

        ## Write integer fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'keepmetric', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'restart', format = 'Integer')

        ## Write string fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'fieldname', format = 'String')
