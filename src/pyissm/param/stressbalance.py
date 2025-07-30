import numpy as np
from . import param_utils
from . import class_registry
from .. import execute

@class_registry.register_class
class stressbalance(class_registry.manage_state):
    """
    Stress balance solution parameters class for ISSM.

    This class encapsulates parameters for configuring stress balance simulations in the ISSM (Ice Sheet System Model) framework.
    It controls velocity constraints, convergence criteria, numerical methods, and other parameters for solving 
    the momentum balance equations in ice sheet dynamics.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    spcvx : ndarray, default=nan
        Velocity constraints in x-direction (NaN means no constraint) [m/yr].
    spcvy : ndarray, default=nan
        Velocity constraints in y-direction (NaN means no constraint) [m/yr].
    spcvx_base : ndarray, default=nan
        Base velocity constraints in x-direction (NaN means no constraint) [m/yr].
    spcvy_base : ndarray, default=nan
        Base velocity constraints in y-direction (NaN means no constraint) [m/yr].
    spcvx_shear : ndarray, default=nan
        Shear velocity constraints in x-direction (NaN means no constraint) [m/yr].
    spcvy_shear : ndarray, default=nan
        Shear velocity constraints in y-direction (NaN means no constraint) [m/yr].
    spcvz : ndarray, default=nan
        Velocity constraints in z-direction (NaN means no constraint) [m/yr].
    restol : float, default=1e-4
        Mechanical equilibrium residual convergence criterion.
    reltol : float, default=0.01
        Velocity relative convergence criterion (NaN: not applied).
    abstol : float, default=10
        Velocity absolute convergence criterion (NaN: not applied).
    ishydrologylayer : int, default=0
        Is hydrology layer enabled.
    isnewton : int, default=0
        Numerical method: 0=Picard's fixed point, 1=Newton's method, 2=hybrid.
    FSreconditioning : float, default=1e13
        Full-Stokes reconditioning parameter.
    maxiter : int, default=100
        Maximum number of nonlinear iterations.
    shelf_dampening : float, default=0
        Shelf dampening parameter.
    vertex_pairing : float, default=nan
        Vertex pairing parameter.
    penalty_factor : float, default=3
        Penalty factor for constraint enforcement.
    rift_penalty_lock : float, default=10
        Rift penalty lock parameter.
    rift_penalty_threshold : float, default=0
        Rift penalty threshold parameter.
    referential : float, default=nan
        Referential parameter.
    loadingforce : float, default=nan
        Loading force parameter.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the stressbalance parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the stressbalance parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.stressbalance = pyissm.param.stressbalance()
    md.stressbalance.restol = 1e-5
    md.stressbalance.isnewton = 1
    md.stressbalance.maxiter = 200
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.spcvx = np.nan
        self.spcvy = np.nan
        self.spcvx_base = np.nan
        self.spcvy_base = np.nan
        self.spcvx_shear = np.nan
        self.spcvy_shear = np.nan
        self.spcvz = np.nan
        self.restol = pow(10, -4)
        self.reltol = 0.01
        self.abstol = 10
        self.ishydrologylayer = 0
        self.isnewton = 0
        self.FSreconditioning = pow(10, 13)
        #self.icefront = np.nan -- no longer in use
        self.maxiter = 100
        self.shelf_dampening = 0
        self.vertex_pairing = np.nan
        self.penalty_factor = 3
        self.rift_penalty_lock = 10
        self.rift_penalty_threshold = 0
        self.referential = np.nan
        self.loadingforce = np.nan
        self.requested_outputs = 'List of requested outputs' # Default = ['default'] (Vx, Vy, Vel, Pressure)

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   StressBalance solution parameters:\n'

        s += '      Convergence criteria:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'restol', 'mechanical equilibrium residual convergence criterion'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'reltol', 'velocity relative convergence criterion, NaN: not applied'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'abstol', 'velocity absolute convergence criterion, NaN: not applied'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isnewton', '0: Picard\'s fixed point, 1: Newton\'s method, 2: hybrid'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'maxiter', 'maximum number of nonlinear iterations'))
        s += '\n'
        s += '      boundary conditions:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'spcvx', 'x-axis velocity constraint (NaN means no constraint) [m / yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'spcvy', 'y-axis velocity constraint (NaN means no constraint) [m / yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'spcvz', 'z-axis velocity constraint (NaN means no constraint) [m / yr]'))
        s += '\n'
        s += '      MOLHO boundary conditions:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'spcvx_base', 'x-axis basal velocity constraint (NaN means no constraint) [m / yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'spcvy_base', 'y-axis basal velocity constraint (NaN means no constraint) [m / yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'spcvx_shear', 'x-axis shear velocity constraint (NaN means no constraint) [m / yr]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'spcvy_shear', 'y-axis shear velocity constraint (NaN means no constraint) [m / yr]'))
        s += '\n'
        s += '      Rift options:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rift_penalty_threshold', 'threshold for instability of mechanical constraints'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'rift_penalty_lock', 'number of iterations before rift penalties are locked'))
        s += '\n'
        s += '      Penalty options:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'penalty_factor', 'offset used by penalties: penalty = Kmax * 10^offset'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vertex_pairing', 'pairs of vertices that are penalized'))
        s += '\n'
        s += '      Hydrology layer:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'ishydrologylayer', '(SSA only) 0: no subglacial hydrology layer in driving stress, 1: hydrology layer in driving stress'));
        s += '\n'
        s += '      Other:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'shelf_dampening', 'use dampening for floating ice ? Only for FS model'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'FSreconditioning', 'multiplier for incompressibility equation. Only for FS model'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'referential', 'local referential'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'loadingforce', 'loading force applied on each point [N/m^3]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - stressbalance Class'
        return s


    # Marshall method for saving the stressbalance parameters
    def marshall_class(self, prefix, md, fid):
        """
        Marshall the stressbalance parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.

        Returns
        -------
        None
        """

        # White DoubleMat fields (mattype = 1)
        fieldnames = ['spcvx', 'spcvy', 'spcvz']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

        # Write Double fields
        fieldnames = ['restol', 'reltol', 'FSreconditioning', 'penalty_factor']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Double')

        # Write Integer fields
        fieldnames = ['isnewton', 'maxiter', 'shelf_dampening', 'rift_penalty_lock', 'rift_penalty_threshold']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Integer')

        # Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'abstol', format = 'Double', scale = 1. / md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'ishydrologylayer', format = 'Boolean')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'referential', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'vertex_pairing', format = 'DoubleMat', mattype = 3)

        ## Write conditional fields (loading force and MOLHO)
        # Loading force
        if isinstance(self.loadingforce, (list, tuple, np.ndarray)) and np.size(self.loadingforce, 1) == 3:
            execute.WriteData(fid, prefix, name = 'md.stressbalance.loadingforcex', data = self.loadingforce[:, 0], format = 'DoubleMat', mattype = 1)
            execute.WriteData(fid, prefix, name = 'md.stressbalance.loadingforcey', data = self.loadingforce[:, 1], format = 'DoubleMat', mattype = 1)
            execute.WriteData(fid, prefix, name = 'md.stressbalance.loadingforcez', data = self.loadingforce[:, 2], format = 'DoubleMat', mattype = 1)
        
        # MOLHO
        if md.flowequation.isMOLHO:
            fieldnames = ['spcvx_base', 'spcvy_base', 'spcvx_shear', 'spcvy_shear']
            for field in fieldnames:
                execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'DoubleMat', mattype = 1, scale = 1. / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

        ## TODO: Write requested outputs
