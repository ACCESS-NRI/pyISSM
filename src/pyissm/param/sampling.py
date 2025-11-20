import numpy as np
from . import param_utils
from . import class_registry
from .. import execute

@class_registry.register_class
class sampling(class_registry.manage_state):
    """
    Sampling parameters class for ISSM.

    This class encapsulates parameters for stochastic sampling in the ISSM (Ice Sheet System Model) framework.
    It configures parameters for generating random fields using PDE-based operators and autoregressive processes,
    useful for uncertainty quantification and stochastic forcing applications.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    kappa : float, default=nan
        Coefficient of the identity operator in PDE operator (kappa^2 I - Laplacian)^(alpha/2)(tau).
    tau : float, default=0
        Scaling coefficient of the solution.
    beta : float, default=nan
        Coefficient in Robin boundary conditions (to be defined for robin = 1).
    phi : float, default=nan
        Temporal correlation factor for first-order autoregressive process X_t = phi X_{t-1} + noise (|phi|<1 for stationary process, phi = 1 for random walk process).
    alpha : float, default=2
        Exponent in PDE operator (default: 2.0, BiLaplacian covariance operator).
    robin : int, default=0
        Apply Robin boundary conditions (1 if applied and 0 for homogeneous Neumann boundary conditions).
    seed : int, default=-1
        Seed for pseudorandom number generator (given seed if >=0 and random seed if <0).
    requested_outputs : list, default=[]
        Additional outputs requested (not implemented yet).

    Methods
    -------
    __init__(self, other=None)
        Initializes the sampling parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the sampling parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.sampling = pyissm.param.sampling()
    md.sampling.kappa = 0.1
    md.sampling.alpha = 2.0
    md.sampling.phi = 0.9
    md.sampling.robin = 1
    md.sampling.beta = 0.5
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.kappa = np.nan
        self.tau = 0
        self.beta = np.nan
        self.phi = np.nan
        self.alpha = 2
        self.robin = 0
        self.seed = -1
        self.requested_outputs = []

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Sampling parameters:\n'

        s += '      Parameters of PDE operator (kappa^2 I-Laplacian)^(alpha/2)(tau):\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'kappa', 'coefficient of the identity operator'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'tau', 'scaling coefficient of the solution (default: 1.0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'alpha', 'exponent in PDE operator, (default: 2.0, BiLaplacian covariance operator)'))
        s += '\n'
        s += '      Parameters of Robin boundary conditions nabla () \cdot normvec + beta ():\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'robin', 'Apply Robin boundary conditions (1 if applied and 0 for homogenous Neumann boundary conditions) (default: 0)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'beta', 'Coefficient in Robin boundary conditions (to be defined for robin = 1)'))
        s += '\n'
        s += '      Parameters for first-order autoregressive process (X_t = phi X_{t-1} + noise) (if transient):\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'phi', 'Temporal correlation factor (|phi|<1 for stationary process, phi = 1 for random walk process) (default 0)'))
        s += '\n'
        s += '      Other parameters of stochastic sampler:\n'
        s += '{}\n'.format(param_utils.fielddisplay(self, 'seed', 'Seed for pseudorandom number generator (given seed if >=0 and random seed if <0) (default: -1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested (not implemented yet)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - sampling Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        # Early return if required analysis/solutions are not present
        if 'SamplingAnalysis' not in analyses:
            return md

        param_utils.check_field(md, fieldname = 'sampling.kappa', size = (md.mesh.numberofvertices, ), gt = 0, allow_nan = False, allow_inf = False)
        param_utils.check_field(md, fieldname = 'sampling.tau', gt = 0, numel = 1, allow_nan = False, allow_inf = False)
        param_utils.check_field(md, fieldname = 'sampling.robin', numel = 1, values = [0, 1])
        if md.sampling.robin:
            param_utils.check_field(md, fieldname = 'sampling.beta', size = (md.mesh.numberofvertices, ), gt = 0, allow_nan = False, allow_inf = False)
    
        param_utils.check_field(md, fieldname = 'sampling.alpha', numel = 1, gt = 0, allow_nan = False, allow_inf = False)
        param_utils.check_field(md, fieldname = 'sampling.seed', numel = 1, allow_nan = False, allow_inf = False)
        param_utils.check_field(md, fieldname = 'sampling.requested_outputs', string_list = True)

        return md
    
    # Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self,
                        md = None,
                        return_default_outputs = False):
        """
        Process requested outputs, expanding 'default' to appropriate outputs.

        Parameters
        ----------
        md : ISSM model object, optional
            Model object containing mesh information.
        return_default_outputs : bool, default=False
            Whether to also return the list of default outputs.
            
        Returns
        -------
        outputs : list
            List of output strings with 'default' expanded to actual output names.
        default_outputs : list, optional
            Returned only if `return_default_outputs=True`.
        """

        outputs = []

        ## Set default_outputs
        default_outputs = ['']

        ## Loop through all requested outputs
        for item in self.requested_outputs:
            
            ## Process default outputs
            if item == 'default':
                    outputs.extend(default_outputs)

            ## Append other requested outputs (not defaults)
            else:
                outputs.append(item)

        if return_default_outputs:
            return outputs, default_outputs
        return outputs

    # Marshall method for saving the sampling parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [sampling] parameters to a binary file.

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

        ## Write DoubleMat fields (all consistent formats)
        fieldnames = ['kappa', 'tau', 'beta', 'phi']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'DoubleMat', mattype = 1)

        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'alpha', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'robin', format = 'Boolean')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'seed', format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.sampling.requested_outputs', data = self.process_outputs(md), format = 'StringArray')