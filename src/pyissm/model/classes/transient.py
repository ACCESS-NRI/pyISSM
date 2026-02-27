from pyissm.model.classes import class_utils
from pyissm.model.classes import materials
from pyissm.model.classes import class_registry
from pyissm.model import execute

@class_registry.register_class
class transient(class_registry.manage_state):
    """
    Transient solution parameters class for ISSM.

    This class encapsulates parameters for configuring transient (time-dependent) simulations in the ISSM (Ice Sheet System Model) framework.
    It allows users to enable or disable various physics components and models that can be included in transient simulations,
    such as age tracking, surface mass balance, thermal evolution, and grounding line migration.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    isage : int, default=0
        Indicates if age model is requested in the transient.
    issmb : int, default=1
        Indicates if a surface mass balance solution is used in the transient.
    ismasstransport : int, default=1
        Indicates if a masstransport solution is used in the transient.
    ismmemasstransport : int, default=0
        Indicates whether an MME masstransport solution is used in the transient.
    isoceantransport : int, default=0
        Indicates whether an ocean masstransport solution is used in the transient.
    isstressbalance : int, default=1
        Indicates if a stressbalance solution is used in the transient.
    isthermal : int, default=1
        Indicates if a thermal solution is used in the transient.
    isgroundingline : int, default=0
        Indicates if a groundingline migration is used in the transient.
    isesa : int, default=0
        Indicates whether an elastic adjustment model is used in the transient.
    isdamageevolution : int, default=0
        Indicates whether damage evolution is used in the transient.
    ismovingfront : int, default=0
        Indicates whether a moving front capability is used in the transient.
    ishydrology : int, default=0
        Indicates whether an hydrology model is used.
    isdebris : int, default=0
        Indicates whether a debris model is used.
    issampling : int, default=0
        Indicates whether sampling is used in the transient.
    isslc : int, default=0
        Indicates if a sea level change solution is used in the transient.
    amr_frequency : int, default=0
        Frequency at which mesh is refined in simulations with multiple time_steps.
    isoceancoupling : int, default=0
        Indicates whether coupling with an ocean model is used in the transient (1 for cartesian coordinates, 2 for lat/long coordinates).
    requested_outputs : list, default=['default']
        List of additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the transient parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the transient parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.transient = pyissm.model.classes.transient()
    md.transient.isage = 1
    md.transient.isgroundingline = 1
    md.transient.requested_outputs = ['IceVolume', 'IceVolumeAboveFloatation']
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.isage = 0
        self.issmb = 1
        self.ismasstransport = 1
        self.ismmemasstransport = 0
        self.isoceantransport = 0
        self.isstressbalance = 1
        self.isthermal = 1
        self.isgroundingline = 0
        self.isesa = 0
        self.isdamageevolution = 0
        self.ismovingfront = 0
        self.ishydrology = 0
        self.isdebris = 0
        self.issampling = 0
        self.isslc = 0
        self.amr_frequency = 0
        self.isoceancoupling = 0
        self.requested_outputs = ['default']

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   transient solution parameters:\n'

        s += '{}\n'.format(class_utils._field_display(self, 'isage', 'indicates if age model is requested in the transient'))
        s += '{}\n'.format(class_utils._field_display(self, 'issmb', 'indicates if a surface mass balance solution is used in the transient'))
        s += '{}\n'.format(class_utils._field_display(self, 'ismasstransport', 'indicates if a masstransport solution is used in the transient'))
        s += '{}\n'.format(class_utils._field_display(self, 'ismmemasstransport', 'indicates whether an MME masstransport solution is used in the transient'))
        s += '{}\n'.format(class_utils._field_display(self, 'isoceantransport', 'indicates whether an ocean masstransport solution is used in the transient'))
        s += '{}\n'.format(class_utils._field_display(self, 'isstressbalance', 'indicates if a stressbalance solution is used in the transient'))
        s += '{}\n'.format(class_utils._field_display(self, 'isthermal', 'indicates if a thermal solution is used in the transient'))
        s += '{}\n'.format(class_utils._field_display(self, 'isgroundingline', 'indicates if a groundingline migration is used in the transient'))
        s += '{}\n'.format(class_utils._field_display(self, 'isesa', 'indicates whether an elastic adjustment model is used in the transient'))
        s += '{}\n'.format(class_utils._field_display(self, 'isdamageevolution', 'indicates whether damage evolution is used in the transient'))
        s += '{}\n'.format(class_utils._field_display(self, 'ismovingfront', 'indicates whether a moving front capability is used in the transient'))
        s += '{}\n'.format(class_utils._field_display(self, 'ishydrology', 'indicates whether an hydrology model is used'))
        s += '{}\n'.format(class_utils._field_display(self, 'isdebris', 'indicates whether a debris model is used'))
        s += '{}\n'.format(class_utils._field_display(self, 'issampling', 'indicates whether sampling is used in the transient'))
        s += '{}\n'.format(class_utils._field_display(self, 'isslc', 'indicates if a sea level change solution is used in the transient'))
        s += '{}\n'.format(class_utils._field_display(self, 'isoceancoupling', 'indicates whether coupling with an ocean model is used in the transient (1 for cartesian coordinates, 2 for lat/long coordinates'))
        s += '{}\n'.format(class_utils._field_display(self, 'amr_frequency', 'frequency at which mesh is refined in simulations with multiple time_steps'))
        s += '{}\n'.format(class_utils._field_display(self, 'requested_outputs', 'list of additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - transient Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        #Early return if not a transient solution
        if not solution == 'TransientSolution':
            return md

        class_utils.check_field(md, fieldname = 'transient.isage', scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = 'transient.issmb', scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = 'transient.ismasstransport', scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = 'transient.ismmemasstransport', scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = 'transient.isoceantransport', scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = 'transient.isstressbalance', scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = 'transient.isthermal', scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = 'transient.isgroundingline', scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = 'transient.isesa', scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = 'transient.isdamageevolution', scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = 'transient.ishydrology', scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = 'transient.isdebris', scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = 'transient.issampling', scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = 'transient.ismovingfront', scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = 'transient.isslc', scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = 'transient.isoceancoupling', scalar = True, values = [0, 1, 2])
        class_utils.check_field(md, fieldname = 'transient.amr_frequency', scalar = True, ge = 0, allow_nan = False, allow_inf = False)

        if solution != 'TransientSolution' and md.transient.iscoupling:
            md.check_message("Coupling with ocean can only be done in transient simulations!")
        
        if md.transient.isdamageevolution and not isinstance(md.materials, materials.damageice):
            md.check_message("requesting damage evolution but md.materials is not of class damageice")
        
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
        default_outputs = []

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
    
    def deactivate_all(self):
        """
        Deactivate all transient model components.
        """

        self.isage = 0
        self.issmb = 0
        self.ismasstransport = 0
        self.ismmemasstransport = 0
        self.isoceantransport = 0
        self.isstressbalance = 0
        self.isthermal = 0
        self.isgroundingline = 0
        self.isesa = 0
        self.isdamageevolution = 0
        self.ismovingfront = 0
        self.ishydrology = 0
        self.isdebris = 0
        self.issampling = 0
        self.isslc = 0
        self.requested_outputs = []

        return self

    # Marshall method for saving the transient parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [transient] parameters to a binary file.

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

        ## Write Boolean fields
        fieldnames = ['isage', 'issmb', 'ismasstransport', 'ismmemasstransport',
            'isoceantransport', 'isstressbalance', 'isthermal',
            'isgroundingline', 'isesa', 'isdamageevolution',
            'ismovingfront', 'ishydrology', 'isdebris',
            'issampling', 'isslc']
        for fieldname in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = fieldname, format = 'Boolean')
        
        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'isoceancoupling', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'amr_frequency', format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.transient.requested_outputs', data = self.process_outputs(md), format = 'StringArray')