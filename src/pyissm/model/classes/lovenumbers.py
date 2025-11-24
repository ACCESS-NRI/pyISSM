import numpy as np
from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm.model import execute

@class_registry.register_class
class lovenumbers(class_registry.manage_state):
    """
    Love numbers parameters class for ISSM.

    This class encapsulates Love numbers parameters for the ISSM (Ice Sheet System Model) framework.
    Love numbers describe the elastic and viscoelastic response of the solid Earth to surface loading,
    including ice mass changes. They are essential for modeling sea level change and solid Earth deformation.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    h : float, default=nan
        Load Love number for radial displacement (provided by PREM model).
    k : float, default=nan
        Load Love number for gravitational potential perturbation (provided by PREM model).
    l : float, default=nan
        Load Love number for horizontal displacements (provided by PREM model).
    th : float, default=nan
        Tidal load Love number (degree 2) for radial displacement.
    tk : float, default=nan
        Tidal load Love number (degree 2) for gravitational potential.
    tl : float, default=nan
        Tidal load Love number (degree 2) for horizontal displacement.
    tk2secular : float, default=0
        Secular fluid Love number (degree 2).
    pmtf_colinear : float, default=nan
        Colinear component of the Polar Motion Transfer Function (e.g. x-motion due to x-component perturbation of the inertia tensor).
    pmtf_ortho : float, default=nan
        Orthogonal component of the Polar Motion Transfer Function (couples x and y components, only used for Chandler Wobble).
    timefreq : float, default=nan
        Time/frequency vector [yr or 1/yr].
    istime : int, default=1
        Time (1, default) or frequency love numbers (0).

    Methods
    -------
    __init__(self, other=None)
        Initializes the lovenumbers parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the lovenumbers parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Notes
    -----
    This functionality is not yet fully implemented in the current version.

    Examples
    --------
    md.lovenumbers = pyissm.model.classes.lovenumbers()
    md.lovenumbers.h = h_love_numbers
    md.lovenumbers.k = k_love_numbers
    md.lovenumbers.l = l_love_numbers
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        # Loading love numbers
        self.h = np.nan # Provided by PREM model
        self.k = np.nan # idem
        self.l = np.nan # idem

        # Tidal love numbers for computing rotational feedback
        self.th = np.nan
        self.tk = np.nan
        self.tl = np.nan
        self.tk2secular = 0 # deg 2 secular number
        self.pmtf_colinear = np.nan
        self.pmtf_ortho = np.nan

        # Time/frequency for visco-elastic love numbers
        self.timefreq = np.nan
        self.istime = 1

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '---------------------------------------\n'
        s += '****      NOT YET IMPLEMENTED      ****\n'
        s += '---------------------------------------\n\n'
        s += '   lovenumbers parameters:\n'
        s += '{}\n'.format(class_utils.fielddisplay(self, 'h', 'load Love number for radial displacement'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'k', 'load Love number for gravitational potential perturbation'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'l', 'load Love number for horizontal displacements'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'th', 'tidal load Love number (deg 2)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'tk', 'tidal load Love number (deg 2)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'tl', 'tidal load Love number (deg 2)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'tk2secular', 'secular fluid Love number'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'pmtf_colinear', 'Colinear component of the Polar Motion Transfer Function (e.g. x-motion due to x-component perturbation of the inertia tensor)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'pmtf_ortho', 'Orthogonal component of the Polar Motion Transfer Function (couples x and y components, only used for Chandler Wobble)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'istime', 'time (default: 1) or frequency love numbers (0)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'timefreq', 'time/frequency vector (yr or 1/yr)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'pmtf_colinear', 'Colinear component of the Polar Motion Transfer Function (e.g. x-motion due to x-component perturbation of the inertia tensor)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'pmtf_ortho', 'Orthogonal component of the Polar Motion Transfer Function (couples x and y components, only used for Chandler Wobble)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - lovenumbers Class'
        return s
    
    def check_consistency(self, md, solution, analyses):
        # Early return if required analyses and solution not requested
        if ('SealevelchangeAnalysis' not in analyses) or (solution == 'TransientSolution' and not md.transient.isslc):
            return
        
        class_utils.check_field(md, fieldname = 'solidearth.lovenumbers.h', allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'solidearth.lovenumbers.k', allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'solidearth.lovenumbers.l', allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'solidearth.lovenumbers.th', allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'solidearth.lovenumbers.tk', allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'solidearth.lovenumbers.tl', allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'solidearth.lovenumbers.tk2secular', allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'solidearth.lovenumbers.pmtf_colinear', allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'solidearth.lovenumbers.pmtf_ortho', allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'solidearth.lovenumbers.timefreq', allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'solidearth.lovenumbers.istime', values = [0, 1], allow_nan = False, allow_inf = False)

        # Check that love numbers are provided at the same level of accuracy
        if (self.h.shape[0] != self.k.shape[0]) or (self.h.shape[0] != self.l.shape[0]):
            raise ValueError('pyissm.model.classes.lovenumbers.check_consistency: love numbers should be provided at the same level of accuracy')

        ntf = len(self.timefreq)
        if (np.shape(self.h)[1] != ntf or np.shape(self.k)[1] != ntf or np.shape(self.l)[1] != ntf or np.shape(self.th)[1] != ntf or np.shape(self.tk)[1] != ntf or np.shape(self.tl)[1] != ntf or np.shape(self.pmtf_colinear)[1] != ntf or np.shape(self.pmtf_ortho)[1] != ntf):
            raise ValueError('pyissm.model.classes.lovenumbers.check_consistency: love numbers should have as many time/frequency steps as the time/frequency vector')

        if self.istime and self.timefreq[0] != 0:
            raise ValueError('pyissm.model.classes.lovenumbers.check_consistency: temporal love numbers must start with elastic response, i.e. timefreq[0] = 0')
        return md
    
    # Marshall method for saving the lovenumbers parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [lovenumbers] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
            ISSM model object needed in some cases.

        Returns
        -------
        None
        """

        ## Write DoubleMat fields
        fieldname = ['h', 'k', 'l', 'th', 'tk', 'tl', 'pmtf_colinear', 'pmtf_ortho']
        for field in fieldname:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'DoubleMat', mattype = 1)

        ## Write conditional fields
        if (self.istime):
            scale = md.constants.yts
        else:
            scale = 1.0 / md.constants.yts

        execute.WriteData(fid, prefix, obj = self, fieldname = 'timefreq', format = 'DoubleMat', mattype = 1, scale = scale)

        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'tk2secular', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'istime', format = 'Boolean')
