from . import param_utils
from . import class_registry
from .. import execute

@class_registry.register_class
class issmsettings(class_registry.manage_state):
    """
    ISSM general settings class for ISSM.

    This class encapsulates general settings and configuration parameters for the ISSM (Ice Sheet System Model) framework.
    It controls output frequency, I/O behavior, memory management, solver settings, and other system-wide parameters
    that affect the overall behavior of ISSM simulations.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    results_on_nodes : str, default='List of output'
        List of output for which results will be output for all the nodes of each element. Use 'all' for all output on nodes.
    io_gather : int, default=1
        I/O gathering strategy for result outputs.
    lowmem : int, default=0
        Is the memory limited? (0 or 1).
    output_frequency : int, default=1
        Number of time steps between two saves (e.g., 5 means that results are only saved every 5 time steps).
    sb_coupling_frequency : int, default=1
        Frequency at which StressBalance solver is coupled.
    checkpoint_frequency : int, default=0
        Frequency at which the runs are being recorded, allowing for a restart.
    waitonlock : int, default=2147483647
        Maximum number of minutes to wait for batch results, or return 0 (default is 2^31-1).
    solver_residue_threshold : float, default=1e-6
        Throw an error if solver residue exceeds this value (NaN to deactivate).

    Methods
    -------
    __init__(self, other=None)
        Initializes the issmsettings parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the issmsettings parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.settings = pyissm.param.issmsettings()
    md.settings.output_frequency = 5
    md.settings.lowmem = 1
    md.settings.results_on_nodes = ['Vel', 'Thickness']
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.results_on_nodes = 'List of output'
        self.io_gather = 1
        self.lowmem = 0
        self.output_frequency = 1
        self.sb_coupling_frequency = 1
        self.checkpoint_frequency = 0
        self.waitonlock = pow(2, 31) - 1
        self.solver_residue_threshold = 1e-6

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   general issmsettings parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'results_on_nodes', "list of output for which results will be output for all the nodes of each element, Use 'all' for all output on nodes."))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'io_gather', 'I / O gathering strategy for result outputs (default 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'lowmem', 'is the memory limited ? (0 or 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'output_frequency', 'number of time steps between two saves (e.g., 5 means that results are only saved every 5 time steps)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'sb_coupling_frequency', 'frequency at which StressBalance solver is coupled (default 1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'checkpoint_frequency', 'frequency at which the runs are being recorded, allowing for a restart'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'waitonlock', 'maximum number of minutes to wait for batch results, or return 0'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'solver_residue_threshold', 'throw an error if solver residue exceeds this value (NaN to deactivate)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - issmsettings Class'
        return s

    # Marshall method for saving the settings parameters
    def marshall_class(self, prefix, md, fid):
        """
        Marshall the settings parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.

        Returns
        -------
        None
        """
        
        ## Write the settings parameters to the file
        execute.WriteData(fid, prefix, name = 'md.settings.results_on_nodes', data = self.results_on_nodes, format = 'StringArray')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'io_gather', format = 'Boolean')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'lowmem', format = 'Boolean')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'output_frequency', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'sb_coupling_frequency', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'checkpoint_frequency', format = 'Integer')

        if self.waitonlock > 0:
            execute.WriteData(fid, prefix, name = 'md.settings.waitonlock', data = True, format = 'Boolean')
        else:
            execute.WriteData(fid, prefix, name = 'md.settings.waitonlock', data = False, format = 'Boolean')

        execute.WriteData(fid, prefix, obj = self, fieldname = 'solver_residue_threshold', format = 'Double')