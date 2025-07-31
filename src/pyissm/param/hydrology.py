import numpy as np
from . import param_utils
from . import class_registry
from .. import execute

## ------------------------------------------------------
## hydrology.armapw
## ------------------------------------------------------
@class_registry.register_class
class armapw(class_registry.manage_state):
    """
    ARMAPW hydrology parameters class for ISSM.

    This class defines the default parameters for the ARMA piecewise (armapw) hydrology model in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    num_basins : int, default=0
        Number of different basins [unitless].
    num_params : int, default=0
        Number of different parameters in the piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.).
    num_breaks : int, default=0
        Number of different breakpoints in the piecewise-polynomial (separating num_breaks+1 periods).
    polynomialparams : ndarray, default=np.nan
        Coefficients for the polynomial (const, trend, quadratic, etc.), dimensioned by basins, periods, and orders.
    arma_timestep : float, default=0
        Time resolution of the ARMA model [yr].
    ar_order : int, default=0
        Order of the autoregressive model [unitless].
    ma_order : int, default=0
        Order of the moving-average model [unitless].
    arlag_coefs : ndarray, default=np.nan
        Basin-specific vectors of AR lag coefficients [unitless].
    malag_coefs : ndarray, default=np.nan
        Basin-specific vectors of MA lag coefficients [unitless].
    datebreaks : ndarray, default=np.nan
        Dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr].
    basin_id : ndarray, default=np.nan
        Basin number assigned to each element [unitless].
    monthlyfactors : ndarray, default=np.nan
        Monthly multiplicative factor on the subglacial water pressure, specified per basin (size: [num_basins, 12]).

    Methods
    -------
    __init__(self, other=None)
        Initializes the armapw parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the armapw parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.hydrology = pyissm.param.hydrology.armapw()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.num_basins = 0
        self.num_params = 0
        self.num_breaks = 0
        self.polynomialparams = np.nan
        self.arma_timestep = 0
        self.ar_order = 0
        self.ma_order = 0
        self.arlag_coefs = np.nan
        self.malag_coefs = np.nan
        self.datebreaks = np.nan
        self.basin_id = np.nan
        self.monthlyfactors = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   hydrologyarmapw parameters:\n'

        s += 'subglacial water pressure is calculated as Pw=monthlyfactor[month]*(rho_water*g*bed+Pw_arma) where Pw_arma is the perturbation calculated as an ARMA process\n'
        s += 'polynomialparams includes the constant, linear trend, quadratic trend, etc. of the ARMA process\n'
        s += 'arlag_coefs and malag_coefs include the coefficients of the ARMA process\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'num_basins', 'number of different basins [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'basin_id', 'basin number assigned to each element [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'num_breaks', 'number of different breakpoints in the piecewise-polynomial (separating num_breaks+1 periods)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'num_params', 'number of different parameters in the piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'monthlyfactors', 'monthly multiplicative factor on the subglacial water pressure, specified per basin (size:[num_basins,12])'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'polynomialparams', 'coefficients for the polynomial (const,trend,quadratic,etc.),dim1 for basins,dim2 for periods,dim3 for orders, ex: polyparams=cat(num_params,intercepts,trendlinearcoefs,trendquadraticcoefs)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'datebreaks', 'dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ar_order', 'order of the autoregressive model [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ma_order', 'order of the moving-average model [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'arma_timestep', 'time resolution of the ARMA model [yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'arlag_coefs', 'basin-specific vectors of AR lag coefficients [unitless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'malag_coefs', 'basin-specific vectors of MA lag coefficients [unitless]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - hydrology.armapw Class'
        return s

## ------------------------------------------------------
## hydrology.dc
## ------------------------------------------------------
@class_registry.register_class
class dc(class_registry.manage_state):
    """
    Dual Porous Continuum Equivalent (DC) hydrology parameters class for ISSM.

    This class defines the default parameters for the dual continuum (dc) hydrology model in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    water_compressibility : float, default=5.04e-10
        Compressibility of water [Pa^-1].
    isefficientlayer : int, default=1
        Use efficient drainage system [1: true, 0: false].
    penalty_factor : int, default=3
        Exponent used in the penalisation method [dimensionless].
    penalty_lock : int, default=0
        Stabilize unstable constraints (default 0: no stabilization).
    rel_tol : float, default=1.0e-06
        Tolerance for nonlinear iteration between layers [dimensionless].
    max_iter : int, default=100
        Maximum number of nonlinear iterations.
    steps_per_step : int, default=1
        Number of hydrology steps per time step.
    step_adapt : int, default=0
        Adaptive sub-stepping [1: true, 0: false].
    averaging : int, default=0
        Averaging method for steps (0: Arithmetic, 1: Geometric, 2: Harmonic).
    sedimentlimit_flag : int, default=0
        Type of upper limit for the inefficient layer (0: none, 1: user, 2: hydrostatic, 3: normal stress).
    sedimentlimit : float, default=0
        User-defined upper limit for the inefficient layer [m].
    transfer_flag : int, default=1
        Transfer method between layers (0: none, 1: constant leakage).
    unconfined_flag : int, default=0
        Use unconfined scheme (0: confined only, 1: confined-unconfined).
    leakage_factor : float, default=1.0e-10
        User-defined leakage factor [m].
    basal_moulin_input : ndarray, default=np.nan
        Water flux at a given point [m3 s^-1].
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.
    spcsediment_head : ndarray, default=np.nan
        Sediment water head constraints [m above MSL].
    mask_thawed_node : ndarray, default=np.nan
        Mask for thawed nodes (0: frozen).
    sediment_transmitivity : float, default=8.0e-04
        Sediment transmissivity [m^2/s].
    sediment_compressibility : float, default=1.0e-08
        Sediment compressibility [Pa^-1].
    sediment_porosity : float, default=0.4
        Sediment porosity [dimensionless].
    sediment_thickness : float, default=20.0
        Sediment thickness [m].
    spcepl_head : ndarray, default=np.nan
        EPL water head constraints [m above MSL].
    mask_eplactive_node : ndarray, default=np.nan
        Mask for active EPL nodes (1: active, 0: inactive).
    epl_compressibility : float, default=1.0e-08
        EPL compressibility [Pa^-1].
    epl_porosity : float, default=0.4
        EPL porosity [dimensionless].
    epl_initial_thickness : float, default=1.0
        EPL initial thickness [m].
    epl_thick_comp : int, default=1
        EPL thickness computation flag.
    epl_max_thickness : float, default=5.0
        EPL maximal thickness [m].
    epl_conductivity : float, default=8.0e-02
        EPL conductivity [m^2/s].
    epl_colapse_thickness : float
        EPL collapsing thickness [m] (computed as sediment_transmitivity / epl_conductivity).
    eplflip_lock : int, default=0
        Lock EPL activity to avoid flip-flopping (default 0: no stabilization).

    Methods
    -------
    __init__(self, other=None)
        Initializes the dc parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the dc parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.hydrology = pyissm.param.hydrology.dc()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.water_compressibility = 5.04e-10
        self.isefficientlayer = 1
        self.penalty_factor = 3
        self.penalty_lock = 0
        self.rel_tol = 1.0e-06
        self.max_iter = 100
        self.steps_per_step = 1
        self.step_adapt = 0
        self.averaging = 0
        self.sedimentlimit_flag = 0
        self.sedimentlimit = 0
        self.transfer_flag = 1
        self.unconfined_flag = 0
        self.leakage_factor = 1.0e-10
        self.basal_moulin_input = np.nan
        self.requested_outputs = 'List of requested outputs'
        self.spcsediment_head = np.nan
        self.mask_thawed_node = np.nan
        self.sediment_transmitivity = 8.0e-04
        self.sediment_compressibility = 1.0e-08
        self.sediment_porosity = 0.4
        self.sediment_thickness = 20.0
        self.spcepl_head = np.nan
        self.mask_eplactive_node = np.nan
        self.epl_compressibility = 1.0e-08
        self.epl_porosity = 0.4
        self.epl_initial_thickness = 1.0
        self.epl_thick_comp = 1
        self.epl_max_thickness = 5.0
        self.epl_conductivity = 8.0e-02
        self.epl_colapse_thickness = self.sediment_transmitivity / self.epl_conductivity
        self.eplflip_lock = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   hydrology Dual Porous Continuum Equivalent parameters:\n'

        s += '\t- general parameters\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'water_compressibility', 'compressibility of water [Pa^ - 1]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isefficientlayer', 'do we use an efficient drainage system [1: true 0: false]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'penalty_factor', 'exponent of the value used in the penalisation method [dimensionless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'penalty_lock', 'stabilize unstable constraints that keep zigzagging after n iteration (default is 0, no stabilization)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rel_tol', 'tolerance of the nonlinear iteration for the transfer between layers [dimensionless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'max_iter', 'maximum number of nonlinear iteration'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'steps_per_step', 'number of hydrology steps per time step'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'step_adapt', 'adaptative sub stepping  [1: true 0: false] default is 0'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '{}\n'.format('                   0: Arithmetic (default)')
        s += '{}\n'.format('                   1: Geometric')
        s += '{}\n'.format('                   2: Harmonic')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'basal_moulin_input', 'water flux at a given point [m3 s - 1]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sedimentlimit_flag', 'what kind of upper limit is applied for the inefficient layer'))
        s += '{}\n'.format('                   0: no limit')
        s += '{}\n'.format('                   1: user defined sedimentlimit')
        s += '{}\n'.format('                   2: hydrostatic pressure')
        s += '{}\n'.format('                   3: normal stress')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sedimentlimit', 'user defined upper limit for the inefficient layer [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'transfer_flag', 'what kind of transfer method is applied between the layers'))
        s += '{}\n'.format('                   0: no transfer')
        s += '{}\n'.format('                   1: constant leakage factor: leakage_factor')
        s += '{}\n'.format(param_utils.fielddisplay(self, 'leakage_factor', 'user defined leakage factor [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'unconfined_flag', 'using an unconfined scheme or not (transitory)'))
        s += '{}\n'.format('                   0: Confined only')
        s += '{}\n'.format('                   1: Confined - Unconfined')
        s += '\t- for the sediment layer\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'spcsediment_head', 'sediment water head constraints (NaN means no constraint) [m above MSL]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sediment_compressibility', 'sediment compressibility [Pa^ - 1]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sediment_porosity', 'sediment [dimensionless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sediment_thickness', 'sediment thickness [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sediment_transmitivity', 'sediment transmitivity [m^2 / s]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'mask_thawed_node', 'IDS is deactivaed (0) on frozen nodes'))
        s += '\t- for the epl layer\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'spcepl_head', 'epl water head constraints (NaN means no constraint) [m above MSL]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'mask_eplactive_node', 'active (1) or not (0) EPL'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'epl_compressibility', 'epl compressibility [Pa^ - 1]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'epl_porosity', 'epl [dimensionless]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'epl_max_thickness', 'epl maximal thickness [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'epl_initial_thickness', 'epl initial thickness [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'epl_colapse_thickness', 'epl colapsing thickness [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'epl_thick_comp', 'epl thickness computation flag'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'epl_conductivity', 'epl conductivity [m^2 / s]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'eplflip_lock', 'lock epl activity to avoid flip - floping (default is 0, no stabilization)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - hydrology.dc Class'
        return s

## ------------------------------------------------------
## hydrology.glads
## ------------------------------------------------------
@class_registry.register_class
class glads(class_registry.manage_state):
    """
    GlaDS hydrology parameters class for ISSM.

    This class defines the default parameters for the Glacier Drainage System (GlaDS) hydrology model in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    pressure_melt_coefficient : float, default=7.5e-8
        Pressure melt coefficient (c_t) [K Pa^-1].
    sheet_conductivity : float or ndarray, default=np.nan
        Sheet conductivity (k) [m^(7/4) kg^(-1/2)].
    cavity_spacing : float, default=2.0
        Cavity spacing (l_r) [m].
    bump_height : float or ndarray, default=np.nan
        Typical bump height (h_r) [m].
    omega : float, default=1./2000.
        Transition parameter (omega) [].
    sheet_alpha : float, default=5.0/4.0
        First sheet-flow exponent (alpha_s) [].
    sheet_beta : float, default=3.0/2.0
        Second sheet-flow exponent (beta_s) [].
    rheology_B_base : float or ndarray, default=np.nan
        Ice rheology factor B at base of ice (B) [Pa s^(-1/3)].
    isincludesheetthickness : int, default=0
        Add rho_w*g*h in effective pressure calculation? 1: yes, 0: no.
    creep_open_flag : int, default=1
        Allow cavities to open by creep when N<0? 1: yes, 0: no.

    ischannels : bool, default=False
        Allow for channels? True or False.
    channel_conductivity : float, default=5.e-2
        Channel conductivity (k_c) [m^(3/2) kg^(-1/2)].
    channel_sheet_width : float, default=2.0
        Channel sheet width [m].
    channel_alpha : float, default=5.0/4.0
        First channel-flow exponent (alpha_c) [].
    channel_beta : float, default=3.0/2.0
        Second channel-flow exponent (beta_c) [].

    spcphi : float or ndarray, default=np.nan
        Hydraulic potential Dirichlet constraints [Pa].
    moulin_input : float or ndarray, default=np.nan
        Moulin input (Q_s) [m^3/s].
    neumannflux : float or ndarray, default=np.nan
        Water flux applied along the model boundary [m^2/s].
    englacial_void_ratio : float, default=1.e-5
        Englacial void ratio (e_v).
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.
    melt_flag : int, default=0
        User specified basal melt? 0: no (default), 1: use md.basalforcings.groundedice_melting_rate.
    istransition : int, default=0
        Use standard [0, default] or transition model [1].

    Methods
    -------
    __init__(self, other=None)
        Initializes the glads parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the glads parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.hydrology = pyissm.param.hydrology.glads()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        # Sheet
        self.pressure_melt_coefficient = 7.5e-8
        self.sheet_conductivity = np.nan
        self.cavity_spacing = 2.
        self.bump_height = np.nan
        self.omega = 1./2000.
        self.sheet_alpha = 5.0/4.0
        self.sheet_beta = 3.0/2.0
        self.rheology_B_base = np.nan
        self.isincludesheetthickness = 0
        self.creep_open_flag = 1
        self.rheology_B_base = np.nan

        # Channels
        self.ischannels = False
        self.channel_conductivity = 5.e-2
        self.channel_sheet_width = 2.
        self.channel_alpha = 5.0/4.0
        self.channel_beta = 3.0/2.0

        # Other
        self.spcphi = np.nan
        self.moulin_input = np.nan
        self.neumannflux = np.nan
        self.englacial_void_ratio = 1.e-5
        self.requested_outputs = 'List of requested outputs'
        self.melt_flag = 0
        self.istransition = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   GlaDS (hydrologyglads) solution parameters:\n'

        s += '\t--SHEET\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'pressure_melt_coefficient', 'Pressure melt coefficient (c_t) [K Pa^ - 1]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sheet_conductivity', 'sheet conductivity (k) [m^(7 / 4) kg^(- 1 / 2)]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sheet_alpha', 'First sheet-flow exponent (alpha_s) []'))  # TH
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sheet_beta', 'Second sheet-flow exponent (beta_s) []'))  # TH
        s += '{}\n'.format(param_utils.fielddisplay(self, 'cavity_spacing', 'cavity spacing (l_r) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'bump_height', 'typical bump height (h_r) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'omega', 'transition parameter (omega) []'))  # TH
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rheology_B_base', 'ice rheology factor B at base of ice (B) [Pa s^(-1/3)]'))  # SE
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isincludesheetthickness', 'Do we add rho_w*g*h in effective pressure calculation? 1: yes, 0: no'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'creep_open_flag', 'Do we allow cavities to open by creep when N<0? 1: yes, 0: no'))
        s += '\t--CHANNELS\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ischannels', 'Do we allow for channels? 1: yes, 0: no'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'channel_conductivity', 'channel conductivity (k_c) [m^(3 / 2) kg^(- 1 / 2)]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'channel_sheet_width', 'channel sheet width [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'channel_alpha', 'First channel-flow exponent (alpha_s) []'))  # TH
        s += '{}\n'.format(param_utils.fielddisplay(self, 'channel_beta', 'Second channel-flow exponent (beta_s) []'))  # TH
        s += '\t--OTHER\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'spcphi', 'Hydraulic potential Dirichlet constraints [Pa]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'neumannflux', 'water flux applied along the model boundary (m^2 / s)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'moulin_input', 'moulin input (Q_s) [m^3 / s]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'englacial_void_ratio', 'englacial void ratio (e_v)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'melt_flag', 'User specified basal melt? 0: no (default), 1: use md.basalforcings.groundedice_melting_rate'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'istransition', 'do we use standard [0, default] or transition model [1]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - hydrology.glads Class'
        return s

## ------------------------------------------------------
## hydrology.pism
## ------------------------------------------------------
@class_registry.register_class
class pism(class_registry.manage_state):
    """
    PISM hydrology parameters class for ISSM.

    This class defines the default parameters for the PISM hydrology model in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    drainage_rate : ndarray, default=np.nan
        Fixed drainage rate [mm/yr].
    watercolumn_max : float, default=np.nan
        Maximum water column height [m], recommended default: 2 m.

    Methods
    -------
    __init__(self, other=None)
        Initializes the pism parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the pism parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.hydrology = pyissm.param.hydrology.pism()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.drainage_rate = np.nan
        self.watercolumn_max = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   hydrologypism solution parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'drainage_rate', 'fixed drainage rate [mm / yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'watercolumn_max', 'maximum water column height [m], recommended default: 2 m'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - hydrology.pism Class'
        return s

## ------------------------------------------------------
## hydrology.shakti
## ------------------------------------------------------
@class_registry.register_class
class shakti(class_registry.manage_state):
    """
    Shakti hydrology parameters class for ISSM.

    This class defines the default parameters for the Shakti hydrology model in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    head : float or ndarray, default=np.nan
        Subglacial hydrology water head [m].
    gap_height : float or ndarray, default=np.nan
        Height of gap separating ice from bed [m].
    gap_height_min : float, default=1e-3
        Minimum allowed gap height [m].
    gap_height_max : float, default=1.0
        Maximum allowed gap height [m].
    bump_spacing : float or ndarray, default=np.nan
        Characteristic bedrock bump spacing [m].
    bump_height : float or ndarray, default=np.nan
        Characteristic bedrock bump height [m].
    englacial_input : float or ndarray, default=np.nan
        Liquid water input from englacial to subglacial system [m/yr].
    moulin_input : float or ndarray, default=np.nan
        Liquid water input from moulins (at the vertices) to subglacial system [m^3/s].
    reynolds : float or ndarray, default=np.nan
        Reynolds number.
    spchead : float or ndarray, default=np.nan
        Water head constraints (NaN means no constraint) [m].
    neumannflux : float or ndarray, default=np.nan
        Water flux applied along the model boundary [m^2/s].
    relaxation : float, default=1
        Under-relaxation coefficient for nonlinear iteration.
    storage : float or ndarray, default=np.nan
        Englacial storage coefficient (void ratio).
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the shakti parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the shakti parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.hydrology = pyissm.param.hydrology.shakti()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.head = np.nan
        self.gap_height = np.nan
        self.gap_height_min  = 1e-3
        self.gap_height_max  = 1.
        self.bump_spacing = np.nan
        self.bump_height = np.nan
        self.englacial_input = np.nan
        self.moulin_input = np.nan
        self.reynolds = np.nan
        self.spchead = np.nan
        self.neumannflux = np.nan
        self.relaxation = 1
        self.storage = np.nan
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   hydrologyshakti solution parameters:'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'head', 'subglacial hydrology water head (m)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'gap_height', 'height of gap separating ice to bed (m)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'gap_height_min', 'minimum allowed gap height (m)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'gap_height_max', 'minimum allowed gap height (m)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'bump_spacing', 'characteristic bedrock bump spacing (m)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'bump_height', 'characteristic bedrock bump height (m)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'englacial_input', 'liquid water input from englacial to subglacial system (m / yr)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'moulin_input', 'liquid water input from moulins (at the vertices) to subglacial system (m^3 / s)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'reynolds', 'Reynolds'' number'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'neumannflux', 'water flux applied along the model boundary (m^2 / s)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'spchead', 'water head constraints (NaN means no constraint) (m)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'relaxation', 'under - relaxation coefficient for nonlinear iteration'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'storage', 'englacial storage coefficient (void ratio)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - hydrology.shakti Class'
        return s

## ------------------------------------------------------
## hydrology.shreve
## ------------------------------------------------------
@class_registry.register_class
class shreve(class_registry.manage_state):
    """
    Shreve hydrology parameters class for ISSM.

    This class defines the default parameters for the Shreve hydrology model in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    spcwatercolumn : ndarray, default=np.nan
        Water thickness constraints (NaN means no constraint) [m].
    stabilization : int, default=1
        Artificial diffusivity (default: 1). Can be more than 1 to increase diffusivity.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the shreve parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the shreve parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.hydrology = pyissm.param.hydrology.shreve()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.spcwatercolumn = np.nan
        self.stabilization = 1
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   hydrologyshreve solution parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'spcwatercolumn', 'water thickness constraints (NaN means no constraint) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'stabilization', 'artificial diffusivity (default: 1). can be more than 1 to increase diffusivity.'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - hydrology.shreve Class'
        return s
    
        
    # Marshall method for saving the hydrology parameters
    def marshall_class(self, prefix, md, fid):
        """
        Marshall the hydrology parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.

        Returns
        -------
        None
        """

        ## Write the header
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, name = 'md.hydrology.model', data = 2, format = 'Integer')

        execute.WriteData(fid, prefix, obj = self, fieldname = 'spcwatercolumn', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'stabilization', format = 'Double')

        ## TODO: Implement marshalling logic for requested_outputs
        execute.WriteData(fid, prefix, name = 'md.hydrology.requested_outputs', data = self.requested_outputs, format = 'StringArray')

## ------------------------------------------------------
## hydrology.tws
## ------------------------------------------------------
@class_registry.register_class
class tws(class_registry.manage_state):
    """
    TWS hydrology parameters class for ISSM.

    This class defines the default parameters for the TWS (two water sheet) hydrology model in ISSM.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    spcwatercolumn : ndarray, default=np.nan
        Water thickness constraints (NaN means no constraint) [m].
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the tws parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the tws parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.hydrology = pyissm.param.hydrology.tws()
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.spcwatercolumn = np.nan
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   hydrologytws solution parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'spcwatercolumn', 'water thickness constraints (NaN means no constraint) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - hydrology.tws Class'
        return s