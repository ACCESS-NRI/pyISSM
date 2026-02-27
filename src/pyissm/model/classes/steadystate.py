import numpy as np
from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm.model import execute

@class_registry.register_class
class steadystate(class_registry.manage_state):
    """
    Steady state solution parameters class for ISSM.

    This class encapsulates parameters for configuring steady state simulations in the ISSM (Ice Sheet System Model) framework.
    It allows users to control convergence criteria and iteration limits for finding steady state solutions
    where the ice sheet geometry and flow fields reach equilibrium.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    reltol : float, default=0.01
        Relative tolerance criterion for convergence.
    maxiter : int, default=100
        Maximum number of iterations allowed.
    requested_outputs : list, default=['default']
        Additional requested outputs for the steady state solution.

    Methods
    -------
    __init__(self, other=None)
        Initializes the steadystate parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the steadystate parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.steadystate = pyissm.model.classes.steadystate()
    md.steadystate.reltol = 0.001
    md.steadystate.maxiter = 200
    md.steadystate.requested_outputs = ['IceVolume', 'IceVolumeAboveFloatation']
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.reltol = 0.01
        self.maxiter = 100
        self.requested_outputs = ['default']

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   steadystate solution parameters:\n'

        s += '{}\n'.format(class_utils._field_display(self, 'reltol', 'relative tolerance criterion'))
        s += '{}\n'.format(class_utils._field_display(self, 'maxiter', 'maximum number of iterations'))
        s += '{}\n'.format(class_utils._field_display(self, 'requested_outputs', 'additional requested outputs'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - steadystate Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        # Early return if required analysis/solutions are not present
        if not solution == 'SteadystateSolution':
            return md

        if not md.timestepping.time_step == 0:
            md.check_message("for a steadystate computation, timestepping.time_step must be zero.")

        if np.isnan(md.stressbalance.reltol):
            md.check_message("for a steadystate computation, stressbalance.reltol (relative convergence criterion) must be defined!")

        class_utils.check_field(md, fieldname = 'steadystate.requested_outputs', string_list = True)

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
        _, stressbalance_defaults = md.stressbalance.process_outputs(md, return_default_outputs = True)
        _, thermal_defaults = md.thermal.process_outputs(md, return_default_outputs = True)
        default_outputs = stressbalance_defaults + thermal_defaults

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
    
    # Marshall method for saving the steadystate parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [steadystate] parameters to a binary file.

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

        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'reltol', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'maxiter', format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.steadystate.requested_outputs', data = self.process_outputs(md), format = 'StringArray')


