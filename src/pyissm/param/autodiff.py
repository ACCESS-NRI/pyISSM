import numpy as np
from . import param_utils
from . import class_registry
from .. import execute

@class_registry.register_class
class autodiff(class_registry.manage_state):
    """
    Automatic differentiation parameters class for ISSM.

    This class encapsulates parameters for automatic differentiation (AD) functionality in the ISSM (Ice Sheet System Model) framework.
    It allows users to configure AD settings including dependent and independent variables, memory buffer sizes, 
    and optimization parameters for sensitivity analysis and gradient-based optimization.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    isautodiff : float, default=0.0
        Indicates if automatic differentiation is activated.
    dependents : str, default='List dependents'
        List of dependent variables for AD.
    independents : str, default='List independents'
        List of independent variables for AD.
    driver : str, default='fos_forward'
        ADOLC driver ('fos_forward' or 'fov_forward').
    obufsize : float, default=524288
        Number of operations per buffer (== OBUFSIZE in usrparms.h).
    lbufsize : float, default=524288
        Number of locations per buffer (== LBUFSIZE in usrparms.h).
    cbufsize : float, default=524288
        Number of values per buffer (== CBUFSIZE in usrparms.h).
    tbufsize : float, default=524288
        Number of taylors per buffer (<=TBUFSIZE in usrparms.h).
    gcTriggerMaxSize : float, default=65536
        Free location block sorting/consolidation triggered if allocated locations exceed this value.
    gcTriggerRatio : float, default=2.0
        Free location block sorting/consolidation triggered if the ratio between allocated and used locations exceeds this value.
    tapeAlloc : float, default=15000000
        Iteration count of a priori memory allocation of the AD tape.
    outputTapeMemory : float, default=0.0
        Write AD tape memory statistics to file ad_mem.dat.
    outputTime : float, default=0.0
        Write AD recording and evaluation times to file ad_time.dat.
    enablePreaccumulation : float, default=0.0
        Enable CoDiPack preaccumulation in augmented places.

    Methods
    -------
    __init__(self, other=None)
        Initializes the autodiff parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the autodiff parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.autodiff = pyissm.param.autodiff()
    md.autodiff.isautodiff = 1
    md.autodiff.dependents = ['Vel']
    md.autodiff.independents = ['MaterialsRheologyBbar']
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.isautodiff = 0.
        self.dependents = 'List dependents'
        self.independents = 'List independents'
        self.driver = 'fos_forward'
        self.obufsize = 524288
        self.lbufsize = 524288
        self.cbufsize = 524288
        self.tbufsize = 524288
        self.gcTriggerMaxSize = 65536
        self.gcTriggerRatio = 2.0
        self.tapeAlloc = 15000000
        self.outputTapeMemory = 0.
        self.outputTime = 0.
        self.enablePreaccumulation = 0.

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):  # {{{
        s = '      automatic differentiation parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'isautodiff', "indicates if the automatic differentiation is activated"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'dependents', "list of dependent variables"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'independents', "list of independent variables"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'driver', "ADOLC driver ('fos_forward' or 'fov_forward')"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'obufsize', "Number of operations per buffer (== OBUFSIZE in usrparms.h)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'lbufsize', "Number of locations per buffer (== LBUFSIZE in usrparms.h)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'cbufsize', "Number of values per buffer (== CBUFSIZE in usrparms.h)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'tbufsize', "Number of taylors per buffer (<=TBUFSIZE in usrparms.h)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'gcTriggerRatio', "free location block sorting / consolidation triggered if the ratio between allocated and used locations exceeds gcTriggerRatio"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'gcTriggerMaxSize', "free location block sorting / consolidation triggered if the allocated locations exceed gcTriggerMaxSize)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'tapeAlloc', 'Iteration count of a priori memory allocation of the AD tape'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'outputTapeMemory', 'Write AD tape memory statistics to file ad_mem.dat'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'outputTime', 'Write AD recording and evaluation times to file ad_time.dat'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'enablePreaccumulation', 'Enable CoDiPack preaccumulation in augmented places'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - autodiff Class'
        return s

    # Marshall method for saving the autodiff parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [autodiff] parameters to a binary file.

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

        ## Write control fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'isautodiff', format = 'Boolean')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'driver', format = 'String')

        if self.isautodiff:

            ## Write Double fields
            fieldnames = ['obufsize', 'lbufsize', 'cbufsize', 'tbufsize', 'gcTriggerRatio', 'gcTriggerMaxSize']
            for fieldname in fieldnames:
                execute.WriteData(fid, prefix, obj = self, fieldname = fieldname, format = 'Double')

            ## Write other fields
            execute.WriteData(fid, prefix, obj = self, fieldname = 'tapeAlloc', format = 'Integer')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'outputTapeMemory', format = 'Boolean')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'outputTime', format = 'Boolean')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'enablePreaccumulation', format = 'Boolean')

            ## Write conditional fields
            ## NOTE: Conditional writing taken from $ISSM_DIR/src/m/classes/autodiff.py
            ## 1 - dependent variables
            num_dependent_objects = len(self.dependents)
            execute.WriteData(fid, prefix, name = 'md.autodiff.num_dependent_objects', data = num_dependent_objects, format = 'Integer')

            if num_dependent_objects:
                names = []
                for i, dep in enumerate(self.dependents):
                    names.append(dep.name)

            execute.WriteData(fid, prefix, name = 'md.autodiff.dependent_object_names', data = names, format = 'StringArray')

            ## 2 - independent variables
            num_independent_objects = len(self.independents)
            execute.WriteData(fid, prefix, name = 'md.autodiff.num_independent_objects', data = num_independent_objects, format = 'Integer')

            for indep in self.independents:
                execute.WriteData(fid, prefix, name = 'md.autodiff.independent_name', data = indep.name, format = 'String')
                execute.WriteData(fid, prefix, name = 'md.autodiff.independent_min_parameters', data = indep.min_parameters, format = 'DoubleMat', mattype = 3)
                execute.WriteData(fid, prefix, name = 'md.autodiff.independent_max_parameters', data = indep.max_parameters, format = 'DoubleMat', mattype = 3)
                execute.WriteData(fid, prefix, name = 'md.autodiff.independent_scaling_factor', data = indep.control_scaling_factor, format = 'Double')
                execute.WriteData(fid, prefix, name = 'md.autodiff.independent_control_size', data = indep.control_size, format = 'Integer')

            ## 3 - build index for fos_forward driver
            if self.driver.lower() == 'fos_forward':
                index = 0

            for indep in self.independents:
                if not np.isnan(indep.fos_forward_index):
                    index += indep.fos_forward_index
                    break
                else:
                    if indep.type == 'scalar':
                        index += 1
                    else:
                        index += indep.nods

            index -= 1  # Convert to c-index numbering
            execute.WriteData(fid, prefix, name = 'md.autodiff.fos_forward_index', data = index, format = 'Integer')

            ## 4 - build index for fos_reverse driver
            if self.driver.lower() == 'fos_reverse':
                index = 0

                for dep in self.dependents:
                    if not np.isnan(dep.fos_reverse_index):
                        index += dep.fos_reverse_index
                        break
                    else:
                        index += 1

                index -= 1  # Convert to c-index numbering
                execute.WriteData(fid, prefix, name = 'md.autodiff.fos_reverse_index', data = index, format = 'Integer')


            ## 5 - build index for fov_forward driver
            if self.driver.lower() == 'fov_forward':
                indices = 0

                for indep in self.independents:
                    if indep.fos_forward_index:
                        indices += indep.fov_forward_indices
                        break
                    else:
                        if indep.type == 'scalar':
                            indices += 1
                        else:
                            indices += indep.nods

                index -= 1  # Convert to c-index numbering
                execute.WriteData(fid, prefix, name = 'md.autodiff.fov_forward_indices', data = indices, format = 'IntMat', mattype = 3)

            ## 6 - Deal with mass fluxes
            mass_flux_segments = []
            for dep in self.dependents:
                if dep.name.lower() == 'massflux':
                    mass_flux_segments.append(dep.segments)

            if mass_flux_segments:
                execute.WriteData(fid, prefix, name = 'md.autodiff.mass_flux_segments', data = mass_flux_segments, format = 'MatArray')
                flag = True
            else:
                flag = False
            execute.WriteData(fid, prefix, name = 'md.autodiff.mass_flux_segments_present', data = flag, format = 'Boolean')

            ## Deal with trace keep on
            keep = False

            ## 7 -  From ADOLC userdoc:
            # The optional integer argument keep of trace on determines whether the 
            # numerical values of all active variables are recorded in a buffered 
            # temporary array or file called the taylor stack. This option takes 
            # effect if keep = 1 and prepares the scene for an immediately 
            # following gradient evaluation by a call to a routine implementing the 
            # reverse mode as described in the Section 4 and Section 5.
            #
            if len(self.driver) <= 3:
                keep = False  # there is no "_reverse" string within the driver string
            else:
                if self.driver[3:].lower().startswith("_reverse"):
                    keep = True
                else:
                    keep = False
            execute.WriteData(fid, prefix, name = 'md.autodiff.keep', data = keep, format ='Boolean')

        else:
            execute.WriteData(fid, prefix, name = 'md.autodiff.mass_flux_segments_present', data = False, format = 'Boolean')
            execute.WriteData(fid, prefix, name = 'md.autodiff.keep', data = False, format = 'Boolean')