import numpy as np
from . import param_utils
from . import class_registry

@class_registry.register_class
class damage(class_registry.manage_state):
    """
    Damage mechanics parameters class for ISSM.

    This class encapsulates parameters for damage mechanics in the ISSM (Ice Sheet System Model) framework.
    Damage mechanics models the evolution of cracks and fractures in ice, affecting ice rheology and flow.
    It is particularly important for modeling ice shelf stability, calving, and fracture propagation.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    isdamage : int, default=0
        Is damage mechanics being used? [0 (default) or 1].
    D : float, default=0
        Damage tensor (scalar for now).
    law : float, default=0
        Damage law ['0: analytical', '1: pralong'].
    spcdamage : ndarray, default=nan
        Damage constraints (NaN means no constraint).
    max_damage : float, default=1 - 1e-5
        Maximum possible damage (0 <= max_damage < 1).
    stabilization : float, default=4
        Stabilization method: 0=no stabilization, 1=artificial diffusion, 2=SUPG (not working), 4=flux corrected transport.
    maxiter : float, default=100
        Maximum number of non-linear iterations.
    elementinterp : str, default='P1'
        Interpolation scheme for finite elements ['P1', 'P2'].
    stress_threshold : float, default=1.3e5
        Stress threshold for damage initiation [Pa].
    stress_ubound : float, default=nan
        Stress upper bound for damage healing [Pa].
    kappa : float, default=2.8
        Ductility parameter for stress softening and damage [> 1].
    c1 : float, default=0
        Damage parameter 1.
    c2 : float, default=0
        Damage parameter 2.
    c3 : float, default=0
        Damage parameter 3.
    c4 : float, default=0
        Damage parameter 4.
    healing : float, default=0
        Healing parameter.
    equiv_stress : float, default=0
        Equivalent stress parameter.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the damage parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the damage parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.damage = pyissm.param.damage()
    md.damage.isdamage = 1
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.isdamage = 0
        self.D = 0
        self.law = 0
        self.spcdamage = np.nan
        self.max_damage = 1 - 1e-5
        self.stabilization = 4
        self.maxiter = 100
        self.elementinterp = 'P1'
        self.stress_threshold = 1.3e5
        self.stress_ubound = np.nan
        self.kappa = 2.8
        self.c1 = 0
        self.c2 = 0
        self.c3 = 0
        self.c4 = 0
        self.healing = 0
        self.equiv_stress = 0
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Damage:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isdamage', 'is damage mechanics being used? [0 (default) or 1]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, "D", "damage tensor (scalar for now)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "law", "damage law ['0: analytical', '1: pralong']"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "spcdamage", "damage constraints (NaN means no constraint)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "max_damage", "maximum possible damage (0 <=max_damage < 1)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "stabilization", "0: no stabilization, 1: artificial diffusion, 2: SUPG (not working), 4: flux corrected transport"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "maxiter", "maximum number of non linear iterations"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "elementinterp", "interpolation scheme for finite elements [''P1'', ''P2'']"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "stress_threshold", "stress threshold for damage initiation (Pa)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "stress_ubound", "stress upper bound for damage healing (Pa)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "kappa", "ductility parameter for stress softening and damage [ > 1]"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "c1", "damage parameter 1 "))
        s += '{}\n'.format(param_utils.fielddisplay(self, "c2", "damage parameter 2 "))
        s += '{}\n'.format(param_utils.fielddisplay(self, "c3", "damage parameter 3 "))
        s += '{}\n'.format(param_utils.fielddisplay(self, "c4", "damage parameter 4 "))
        s += '{}\n'.format(param_utils.fielddisplay(self, "healing", "damage healing parameter"))
        s += '{}\n'.format(param_utils.fielddisplay(self, "equiv_stress", "0: von Mises, 1: max principal"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - damage Class'
        return s

