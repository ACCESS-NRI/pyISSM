from . import param_utils
from . import class_registry
from .. import execute

## ------------------------------------------------------
## love.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    """
    Default parameters class for Love number calculations in ISSM.

    This class defines the default parameters for the Love number computation process.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    nfreq : int, default=1
        Number of frequencies sampled (elastic case: 1).
    frequencies : float or array-like, default=0
        Frequencies sampled (0 for elastic case) [Hz].
    sh_nmax : int, default=256
        Maximum spherical harmonic degree.
    sh_nmin : int, default=1
        Minimum spherical harmonic degree.
    g0 : float, default=9.81
        Gravity constant [m/s^2].
    r0 : float, default=6371000.0
        Reference radius [m].
    mu0 : float, default=1e11
        Reference stress [Pa].
    Gravitational_Constant : float, default=6.67259e-11
        Newtonian constant of gravitation [m^3 kg^-1 s^-2].
    chandler_wobble : int, default=0
        Include inertial terms for Chandler wobble in rotational feedback Love numbers.
    allow_layer_deletion : int, default=1
        Allow migration of the integration boundary with increasing spherical harmonic degree.
    underflow_tol : float, default=1e-16
        Threshold for deep to surface Love number ratio to trigger deletion of layers.
    pw_threshold : float, default=1e-3
        Threshold for bypassing the Post-Widder transform for time-dependent Love numbers.
    min_integration_steps : int, default=500
        Minimum number of radial steps per layer.
    max_integration_dr : float, default=1e4
        Maximum length of radial steps [m].
    integration_scheme : int, default=2
        Integration scheme identifier.
    istemporal : int, default=1
        1 for time-dependent Love numbers, 0 for frequency-dependent or elastic.
    n_temporal_iterations : int, default=7
        Number of spectral samples per time step (inverse Laplace transform).
    time : float or array-like, default=0
        Time vector for deformation [s].
    love_kernels : int, default=0
        Compute Love numbers at depth.
    forcing_type : int, default=11
        Nature and depth of the forcing for the Love number calculation.
    inner_core_boundary : int, default=1
        Interface index for inner core boundary forcing.
    core_mantle_boundary : int, default=2
        Interface index for core mantle boundary forcing.
    complex_computation : int, default=0
        Return Love numbers as real (0) or complex (1).
    quad_precision : int, default=0
        Use 32-digit precision for computation and Post-Widder transform.
    debug : int, default=0
        Output yi system matrix prior to solving.
    hypergeom_table1 : int, default=1
        Table 1 for hypergeometric function (EBM rheology).
    hypergeom_table2 : int, default=1
        Table 2 for hypergeometric function (EBM rheology).
    hypergeom_nalpha : int, default=1
        Length of hypergeometric table (EBM rheology).
    hypergeom_nz : int, default=1
        Width of hypergeometric table (EBM rheology).
    hypergeom_z : int, default=0
        Abscissa for hypergeometric table (EBM rheology).

    Methods
    -------
    __init__(self, other=None)
        Initializes the default Love number parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the Love number parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.love = pyissm.param.love.default()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.nfreq = 1
        self.frequencies = 0
        self.sh_nmax = 256
        self.sh_nmin = 1
        self.g0 = 9.81
        self.r0 = 6371 * 1e3
        self.mu0 = 1e11
        self.Gravitational_Constant = 6.67259e-11
        self.chandler_wobble = 0
        self.allow_layer_deletion = 1
        self.underflow_tol = 1e-16
        self.pw_threshold = 1e-3
        self.min_integration_steps = 500
        self.max_integration_dr = 1e4
        self.integration_scheme = 2
        self.istemporal = 1
        self.n_temporal_iterations = 7
        self.time = 0
        self.love_kernels = 0
        self.forcing_type = 11
        self.inner_core_boundary = 1
        self.core_mantle_boundary = 2
        self.complex_computation = 0
        self.quad_precision = 0
        self.debug = 0
        self.hypergeom_table1 = 1
        self.hypergeom_table2 = 1
        self.hypergeom_nalpha = 1
        self.hypergeom_nz = 1
        self.hypergeom_z = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   love parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'nfreq', 'number of frequencies sampled (default: 1, elastic) [Hz]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'frequencies', 'frequencies sampled (convention defaults to 0 for the elastic case) [Hz]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sh_nmax', 'maximum spherical harmonic degree (default: 256, .35 deg, or 40 km at equator)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sh_nmin', 'minimum spherical harmonic degree (default: 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'g0', 'adimensioning constant for gravity (default: 10) [m/s^2]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'r0', 'adimensioning constant for radius (default: 6371*10^3) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'mu0', 'adimensioning constant for stress (default: 10^11) [Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Gravitational_Constant', 'Newtonian constant of gravitation (default: 6.67259e-11 [m^3 kg^-1 s^-2])'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'chandler_wobble', 'includes the inertial terms for the chandler wobble in the rotational feedback love numbers, only for forcing_type=11 (default: 0) (/!\\ 1 has not been validated yet)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'allow_layer_deletion', 'allow for migration of the integration boundary with increasing spherical harmonics degree (default: 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'underflow_tol', 'threshold of deep to surface love number ratio to trigger the deletion of layers (default: 1e-16)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'pw_threshold', 'if relative variation across frequencies is smaller than this ratio, the post-widder transform for time-dependent love numbers is bypassed (default (1e-3)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'min_integration_steps', 'minimum number of radial steps to propagate the yi system from the bottom to the top of each layer (default: 500)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'max_integration_dr', 'maximum length of radial steps to propagate the yi system from the bottom to the top of each layer (default: 10e3) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'istemporal', ['1 for time-dependent love numbers, 0 for frequency-dependent or elastic love numbers (default: 1)', 'If 1: use fourierlove function build_frequencies_from_time to meet consistency']))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'n_temporal_iterations', 'max number of iterations in the inverse Laplace transform. Also the number of spectral samples per time step requested (default: 7)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'time', 'time vector for deformation if istemporal (default: 0) [s]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'love_kernels', 'compute love numbers at depth? (default: 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'forcing_type', 'integer indicating the nature and depth of the forcing for the Love number calculation (default: 11):'))
        s += '{}\n'.format('                                                     1:  Inner core boundary -- Volumic Potential')
        s += '{}\n'.format('                                                     2:  Inner core boundary -- Pressure')
        s += '{}\n'.format('                                                     3:  Inner core boundary -- Loading')
        s += '{}\n'.format('                                                     4:  Inner core boundary -- Tangential traction')
        s += '{}\n'.format('                                                     5:  Core mantle boundary -- Volumic Potential')
        s += '{}\n'.format('                                                     6:  Core mantle boundary -- Pressure')
        s += '{}\n'.format('                                                     7:  Core mantle boundary -- Loading')
        s += '{}\n'.format('                                                     8:  Core mantle boundary -- Tangential traction')
        s += '{}\n'.format('                                                     9:  Surface -- Volumic Potential')
        s += '{}\n'.format('                                                     10: Surface -- Pressure')
        s += '{}\n'.format('                                                     11: Surface -- Loading')
        s += '{}\n'.format('                                                     12: Surface -- Tangential traction ')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'inner_core_boundary', 'interface index in materials.radius locating forcing. Only used for forcing_type 1--4 (default: 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'core_mantle_boundary', 'interface index in materials.radius locating forcing. Only used for forcing_type 5--8 (default: 2)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'complex_computation', 'return love numbers as 0: real (useful for elastic or temporal forms), 1: complex numbers (useful for Fourier spectral form) (default: 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'quad_precision', 'toggle computation love numbers and post-widder transform with 32 digit precision, useful for temporal form (default: 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'debug', 'outputs yi system matrix prior to solving (default: 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'hypergeom_table1', 'table 1 for hypergeometric function, only for EBM rheology (default: [1])'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'hypergeom_table2', 'table 2 for hypergeometric function, only for EBM rheology (default: [1])'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'hypergeom_nalpha', 'length of hypergeometric table, only for EBM rheology (default: 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'hypergeom_nz', 'width of hypergeometric table, only for EBM rheology (default: 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'hypergeom_z','abscissa for hypergeometric table, only for EBM rheology (default: [0])'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - love Class'
        return s
    
    # Marshall method for saving the love.default parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [love.default] parameters to a binary file.

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

        ## Write Integer fields
        fieldnames = ['nfreq', 'sh_nmax', 'sh_nmin', 'min_integration_steps', 'integration_scheme',
                      'n_temporal_iterations', 'forcing_type', 'inner_core_boundary', 'core_mantle_boundary',
                      'hypergeom_nalpha', 'hypergeom_nz']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj=self, fieldname=field, format='Integer')

        ## Write Double fields
        fieldnames = ['g0', 'r0', 'mu0', 'Gravitational_Constant', 'underflow_tol', 'pw_threshold',
                      'max_integration_dr']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj=self, fieldname=field, format='Double')

        ## Write Boolean fields
        fieldnames = ['chandler_wobble', 'allow_layer_deletion', 'istemporal', 'complex_computation',
                      'quad_precision', 'love_kernels', 'debug']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj=self, fieldname=field, format='Boolean')

        ## Write DoubleMat fields
        execute.WriteData(fid, prefix, obj=self, fieldname='frequencies', format='DoubleMat', mattype = 3)
        execute.WriteData(fid, prefix, obj=self, fieldname='hypergeom_table1', format='DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj=self, fieldname='hypergeom_table2', format='DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj=self, fieldname='hypergeom_z', format='DoubleMat', mattype = 3)

        
## ------------------------------------------------------
## love.fourier
## ------------------------------------------------------
@class_registry.register_class
class fourier(class_registry.manage_state):
    """
    Default parameters class for Fourier Love number calculations in ISSM.

    This class defines the default parameters for the Fourier (frequency-dependent) Love number computation process.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    nfreq : int, default=1
        Number of frequencies sampled (elastic case: 1).
    frequencies : float or array-like, default=0
        Frequencies sampled (0 for elastic case) [Hz].
    sh_nmax : int, default=256
        Maximum spherical harmonic degree.
    sh_nmin : int, default=1
        Minimum spherical harmonic degree.
    g0 : float, default=9.81
        Gravity constant [m/s^2].
    r0 : float, default=6371000.0
        Reference radius [m].
    mu0 : float, default=1e11
        Reference stress [Pa].
    Gravitational_Constant : float, default=6.67259e-11
        Newtonian constant of gravitation [m^3 kg^-1 s^-2].
    chandler_wobble : int, default=0
        Include inertial terms for Chandler wobble in rotational feedback Love numbers.
    allow_layer_deletion : int, default=1
        Allow migration of the integration boundary with increasing spherical harmonic degree.
    underflow_tol : float, default=1e-16
        Threshold for deep to surface Love number ratio to trigger deletion of layers.
    pw_threshold : float, default=1e-3
        Threshold for bypassing the Post-Widder transform for time-dependent Love numbers.
    integration_steps_per_layer : int, default=100
        Number of radial steps per layer.
    istemporal : int, default=0
        1 for time-dependent Love numbers, 0 for frequency-dependent or elastic.
    n_temporal_iterations : int, default=8
        Number of spectral samples per time step (inverse Laplace transform).
    time : float or array-like, default=0
        Time vector for deformation [s].
    love_kernels : int, default=0
        Compute Love numbers at depth.
    forcing_type : int, default=11
        Nature and depth of the forcing for the Love number calculation.
    inner_core_boundary : int, default=1
        Interface index for inner core boundary forcing.
    core_mantle_boundary : int, default=2
        Interface index for core mantle boundary forcing.
    complex_computation : int, default=0
        Return Love numbers as real (0) or complex (1).

    Methods
    -------
    __init__(self, other=None)
        Initializes the Fourier Love number parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the Fourier Love number parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.love = pyissm.param.love.fourier()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.nfreq = 1
        self.frequencies = 0
        self.sh_nmax = 256
        self.sh_nmin = 1
        self.g0 = 9.81
        self.r0 = 6371 * 1e3
        self.mu0 = 1e11
        self.Gravitational_Constant = 6.67259e-11
        self.chandler_wobble = 0
        self.allow_layer_deletion = 1
        self.underflow_tol = 1e-16
        self.pw_threshold = 1e-3
        self.integration_steps_per_layer = 100
        self.istemporal = 0
        self.n_temporal_iterations = 8
        self.time = 0
        self.love_kernels = 0
        self.forcing_type = 11
        self.inner_core_boundary = 1
        self.core_mantle_boundary = 2
        self.complex_computation = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   fourier love parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'nfreq', 'number of frequencies sampled (default: 1, elastic) [Hz]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'frequencies', 'frequencies sampled (convention defaults to 0 for the elastic case) [Hz]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sh_nmax', 'maximum spherical harmonic degree (default: 256, .35 deg, or 40 km at equator)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sh_nmin', 'minimum spherical harmonic degree (default: 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'g0', 'adimensioning constant for gravity (default: 10) [m/s^2]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'r0', 'adimensioning constant for radius (default: 6371*10^3) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'mu0', 'adimensioning constant for stress (default: 10^11) [Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'allow_layer_deletion', 'allow for migration of the integration boundary with increasing spherical harmonics degree (default: 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'Gravitational_Constant', 'Newtonian constant of gravitation (default: 6.67259e-11 [m^3 kg^-1 s^-2])'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'chandler_wobble', 'includes the inertial terms for the chandler wobble in the rotational feedback love numbers, only for forcing_type=11 (default: 0) (/!\\ 1 has not been validated yet)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'allow_layer_deletion', 'allow for migration of the integration boundary with increasing spherical harmonics degree (default: 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'pw_threshold', 'if relative variation across frequencies is smaller than this ratio, the post-widder transform for time-dependent love numbers is bypassed (default (1e-3)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'chandler_wobble', 'includes the inertial terms for the chandler wobble in the rotational feedback love numbers, only for forcing_type=11 (default: 0) (/!\\ 1 is untested)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'integration_steps_per_layer', 'number of radial steps to propagate the yi system from the bottom to the top of each layer (default: 100)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'istemporal', ['1 for time-dependent love numbers, 0 for frequency-dependent or elastic love numbers (default: 0)', 'If 1: use fourierlove function build_frequencies_from_time to meet consistency']))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'n_temporal_iterations', 'max number of iterations in the inverse Laplace transform. Also the number of spectral samples per time step requested (default: 8)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'time', 'time vector for deformation if istemporal (default: 0) [s]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'love_kernels', 'compute love numbers at depth? (default: 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'forcing_type', 'integer indicating the nature and depth of the forcing for the Love number calculation (default: 11):'))
        s += '{}\n'.format('                                                     1:  Inner core boundary -- Volumic Potential')
        s += '{}\n'.format('                                                     2:  Inner core boundary -- Pressure')
        s += '{}\n'.format('                                                     3:  Inner core boundary -- Loading')
        s += '{}\n'.format('                                                     4:  Inner core boundary -- Tangential traction')
        s += '{}\n'.format('                                                     5:  Core mantle boundary -- Volumic Potential')
        s += '{}\n'.format('                                                     6:  Core mantle boundary -- Pressure')
        s += '{}\n'.format('                                                     7:  Core mantle boundary -- Loading')
        s += '{}\n'.format('                                                     8:  Core mantle boundary -- Tangential traction')
        s += '{}\n'.format('                                                     9:  Surface -- Volumic Potential')
        s += '{}\n'.format('                                                     10: Surface -- Pressure')
        s += '{}\n'.format('                                                     11: Surface -- Loading')
        s += '{}\n'.format('                                                     12: Surface -- Tangential traction ')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'inner_core_boundary', 'interface index in materials.radius locating forcing. Only used for forcing_type 1--4 (default: 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'core_mantle_boundary', 'interface index in materials.radius locating forcing. Only used for forcing_type 5--8 (default: 2)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - love.fourier Class'
        return s
    
    # Marshall method for saving the love.fourier parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [love.fourier] parameters to a binary file.

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

        ## Write Integer fields
        fieldnames = ['nfreq', 'sh_nmax', 'sh_nmin', 'integration_steps_per_layer', 'n_temporal_iterations',
                      'forcing_type', 'inner_core_boundary', 'core_mantle_boundary']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj=self, fieldname=field, format='Integer')

        ## Write Double fields
        fieldnames = ['g0', 'r0', 'mu0', 'Gravitational_Constant', 'underflow_tol', 'pw_threshold']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj=self, fieldname=field, format='Double')

        ## Write Boolean fields
        fieldnames = ['chandler_wobble', 'allow_layer_deletion', 'istemporal', 'complex_computation',
                      'love_kernels']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj=self, fieldname=field, format='Boolean')

        ## Write DoubleMat fields
        execute.WriteData(fid, prefix, obj=self, fieldname='frequencies', format='DoubleMat', mattype = 3)