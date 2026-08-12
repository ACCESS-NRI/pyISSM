import functools
from importlib.resources import files

import numpy as np
import pandas as pd

from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm.model import execute

_LOVE_NUMBER_TYPES = [
    'loadingverticaldisplacement',
    'loadinggravitationalpotential',
    'loadinghorizontaldisplacement',
    'tidalverticaldisplacement',
    'tidalgravitationalpotential',
    'tidalhorizontaldisplacement',
]

# CF-frame degree-1 overrides (Blewitt 2003, JGR). Only the three
# non-tidal displacement/potential types have a defined CF correction.
_CF_DEGREE_1_OVERRIDE = {
    'loadingverticaldisplacement': -0.269,
    'loadinggravitationalpotential': 0.021,
    'loadinghorizontaldisplacement': 0.134,
}


@functools.lru_cache(maxsize=1)
def _load_love_table():
    """
    Load the packaged PREM Love-number table (degrees 0..10000, columns
    h, k, l, th, tk, tl) as a (10001, 6) ndarray. Cached after first load.
    """
    path = files('pyissm').joinpath('assets/love_numbers.csv')
    with path.open('r') as fh:
        return pd.read_csv(fh).to_numpy()


def get_love_numbers(love_type, referenceframe='CM', maxdeg=1000):
    """
    Look up a Love-number series from the packaged PREM table.

    Parameters
    ----------
    love_type : str
        One of 'loadingverticaldisplacement', 'loadinggravitationalpotential',
        'loadinghorizontaldisplacement', 'tidalverticaldisplacement',
        'tidalgravitationalpotential', 'tidalhorizontaldisplacement'.
    referenceframe : str, default='CM'
        'CM' (center of mass) or 'CF' (center of figure). CF applies a
        degree-1 override (Blewitt 2003) for the three types that have one.
    maxdeg : int, default=1000
        Maximum spherical-harmonic degree to return; must be <= 10000.

    Returns
    -------
    numpy.ndarray
        1-D array of length maxdeg + 1 (degrees 0..maxdeg).
    """
    if maxdeg > 10000:
        raise ValueError('pyissm.model.classes.lovenumbers.get_love_numbers: '
                          'PREM Love numbers are only tabulated up to degree 10000')
    if maxdeg < 0:
        raise ValueError('pyissm.model.classes.lovenumbers.get_love_numbers: '
                          'maxdeg must be >= 0')
    if love_type not in _LOVE_NUMBER_TYPES:
        raise ValueError('pyissm.model.classes.lovenumbers.get_love_numbers: '
                          f'love_type must be one of {_LOVE_NUMBER_TYPES}')
    if referenceframe not in ('CM', 'CF'):
        raise ValueError("pyissm.model.classes.lovenumbers.get_love_numbers: "
                          "referenceframe must be 'CM' or 'CF'")

    col = _LOVE_NUMBER_TYPES.index(love_type)
    series = _load_love_table()[:maxdeg + 1, col].copy()

    if referenceframe == 'CF' and love_type in _CF_DEGREE_1_OVERRIDE and maxdeg >= 1:
        series[1] = _CF_DEGREE_1_OVERRIDE[love_type]

    return series


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
    maxdeg : int, optional
        If given, populates all Love-number fields (elastic, single-epoch) from the packaged PREM table up to this spherical-harmonic degree (<= 10000). Takes precedence over `other`.
    referenceframe : str, default='CM'
        'CM' (center of mass) or 'CF' (center of figure). Only used when `maxdeg` is given.

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
    __init__(self, other=None, maxdeg=None, referenceframe='CM')
        Initializes the lovenumbers parameters, optionally inheriting from another instance and/or populating from the packaged PREM table.
    __repr__(self)
        Returns a detailed string representation of the lovenumbers parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.lovenumbers = pyissm.model.classes.lovenumbers(maxdeg=10000)

    md.lovenumbers = pyissm.model.classes.lovenumbers()
    md.lovenumbers.h = h_love_numbers
    md.lovenumbers.k = k_love_numbers
    md.lovenumbers.l = l_love_numbers
    """

    # Initialise with default parameters
    def __init__(self, other = None, maxdeg = None, referenceframe = 'CM'):
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

        # Populate from the packaged PREM table, taking precedence over `other`
        if maxdeg is not None:
            self._set_default_parameters(maxdeg, referenceframe)

    def _set_default_parameters(self, maxdeg, referenceframe):
        """
        Populate all Love-number fields from the packaged PREM table up to
        `maxdeg`, for the elastic (single-epoch) case.

        Parameters
        ----------
        maxdeg : int
            Maximum spherical-harmonic degree, <= 10000.
        referenceframe : str
            'CM' or 'CF'.
        """
        self.h = get_love_numbers('loadingverticaldisplacement', referenceframe, maxdeg).reshape(-1, 1)
        self.k = get_love_numbers('loadinggravitationalpotential', referenceframe, maxdeg).reshape(-1, 1)
        self.l = get_love_numbers('loadinghorizontaldisplacement', referenceframe, maxdeg).reshape(-1, 1)
        self.th = get_love_numbers('tidalverticaldisplacement', referenceframe, maxdeg).reshape(-1, 1)
        self.tk = get_love_numbers('tidalgravitationalpotential', referenceframe, maxdeg).reshape(-1, 1)
        self.tl = get_love_numbers('tidalhorizontaldisplacement', referenceframe, maxdeg).reshape(-1, 1)
        self.tk2secular = 0.942
        self.pmtf_colinear = np.array([[0.0]])
        self.pmtf_ortho = np.array([[0.0]])
        if maxdeg >= 2:
            self.pmtf_colinear = ((1.0 + self.k[2, :]) /
                                   (1.0 - self.tk[2, :] / self.tk2secular)).reshape(-1, 1)
        self.istime = 1
        self.timefreq = np.zeros(1)

    # Define repr
    def __repr__(self):
        s = '   lovenumbers parameters:\n'
        s += '{}\n'.format(class_utils._field_display(self, 'h', 'load Love number for radial displacement'))
        s += '{}\n'.format(class_utils._field_display(self, 'k', 'load Love number for gravitational potential perturbation'))
        s += '{}\n'.format(class_utils._field_display(self, 'l', 'load Love number for horizontal displacements'))
        s += '{}\n'.format(class_utils._field_display(self, 'th', 'tidal load Love number (deg 2)'))
        s += '{}\n'.format(class_utils._field_display(self, 'tk', 'tidal load Love number (deg 2)'))
        s += '{}\n'.format(class_utils._field_display(self, 'tl', 'tidal load Love number (deg 2)'))
        s += '{}\n'.format(class_utils._field_display(self, 'tk2secular', 'secular fluid Love number'))
        s += '{}\n'.format(class_utils._field_display(self, 'pmtf_colinear', 'Colinear component of the Polar Motion Transfer Function (e.g. x-motion due to x-component perturbation of the inertia tensor)'))
        s += '{}\n'.format(class_utils._field_display(self, 'pmtf_ortho', 'Orthogonal component of the Polar Motion Transfer Function (couples x and y components, only used for Chandler Wobble)'))
        s += '{}\n'.format(class_utils._field_display(self, 'istime', 'time (default: 1) or frequency love numbers (0)'))
        s += '{}\n'.format(class_utils._field_display(self, 'timefreq', 'time/frequency vector (yr or 1/yr)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - lovenumbers Class'
        return s
    
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [lovenumbers.lovenumbers] parameters.

        Parameters
        ----------
        md : :class:`pyissm.model.Model`
            The model object to check.
        solution : :class:`str`
            The solution name to check.
        analyses : list of :class:`str`
            List of analyses to check consistency for.

        Returns 
        -------
        md : :class:`pyissm.model.Model`
            The model object with any consistency errors noted.
        """
        # Early return if required analyses and solution not requested
        if ('SealevelchangeAnalysis' not in analyses) or (solution == 'TransientSolution' and not md.transient.isslc):
            return
        
        class_utils._check_field(md, fieldname = 'solidearth.lovenumbers.h', allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'solidearth.lovenumbers.k', allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'solidearth.lovenumbers.l', allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'solidearth.lovenumbers.th', allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'solidearth.lovenumbers.tk', allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'solidearth.lovenumbers.tl', allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'solidearth.lovenumbers.tk2secular', allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'solidearth.lovenumbers.pmtf_colinear', allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'solidearth.lovenumbers.pmtf_ortho', allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'solidearth.lovenumbers.timefreq', allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'solidearth.lovenumbers.istime', values = [0, 1], allow_nan = False, allow_inf = False)

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
        Marshall [lovenumbers.lovenumbers] parameters to a binary file.

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
            execute._write_model_field(fid, prefix, obj = self, fieldname = field, format = 'DoubleMat', mattype = 1)

        ## Write conditional fields
        if (self.istime):
            scale = md.constants.yts
        else:
            scale = 1.0 / md.constants.yts

        execute._write_model_field(fid, prefix, obj = self, fieldname = 'timefreq', format = 'DoubleMat', mattype = 1, scale = scale)

        ## Write other fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'tk2secular', format = 'Double')
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'istime', format = 'Boolean')
