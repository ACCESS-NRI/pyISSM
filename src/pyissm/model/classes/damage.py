import numpy as np
from pyissm.model.classes import class_utils, class_registry
from pyissm.model import execute, mesh

@class_registry.register_class
class damage(class_registry.manage_state):
    """
    Damage mechanics class for ISSM.

    This class contains parameters for damage mechanics in the ISSM framework.
    Damage mechanics models the evolution of cracks and fractures in ice, affecting
    ice rheology and flow. It is particularly important for modeling ice shelf stability,
    calving, and fracture propagation.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    isdamage : :class:`int`, default=0
        Is damage mechanics being used? [0 (default) or 1].
    D : :class:`float`, default=0
        Damage tensor (scalar for now).
    law : :class:`float`, default=0
        Damage law ['0: analytical', '1: pralong'].
    spcdamage : :class:`numpy.ndarray`, default=np.nan
        Damage constraints (NaN means no constraint).
    max_damage : :class:`float`, default=1 - 1e-5
        Maximum possible damage (0 <= max_damage < 1).
    stabilization : :class:`float`, default=4
        Stabilization method: 0=no stabilization, 1=artificial diffusion, 2=SUPG (not working), 4=flux corrected transport.
    maxiter : :class:`float`, default=100
        Maximum number of non-linear iterations.
    elementinterp : :class:`str`, default='P1'
        Interpolation scheme for finite elements ['P1', 'P2'].
    stress_threshold : :class:`float`, default=1.3e5
        Stress threshold for damage initiation [Pa].
    stress_ubound : :class:`float`, default=np.nan
        Stress upper bound for damage healing [Pa].
    kappa : :class:`float`, default=2.8
        Ductility parameter for stress softening and damage [> 1].
    c1 : :class:`float`, default=0
        Damage parameter 1.
    c2 : :class:`float`, default=0
        Damage parameter 2.
    c3 : :class:`float`, default=0
        Damage parameter 3.
    c4 : :class:`float`, default=0
        Damage parameter 4.
    healing : :class:`float`, default=0
        Healing parameter.
    equiv_stress : :class:`float`, default=0
        Equivalent stress parameter.
    requested_outputs : :class:`list`, default=['default']
        Additional outputs requested.

    Examples
    --------
    .. code-block:: python
    
        >>> md.damage = pyissm.model.classes.damage()
        >>> md.damage.isdamage = 1
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
        self.requested_outputs = ['default']

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Damage:\n'
        s += '{}\n'.format(class_utils.fielddisplay(self, 'isdamage', 'is damage mechanics being used? [0 (default) or 1]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, "D", "damage tensor (scalar for now)"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "law", "damage law ['0: analytical', '1: pralong']"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "spcdamage", "damage constraints (NaN means no constraint)"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "max_damage", "maximum possible damage (0 <=max_damage < 1)"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "stabilization", "0: no stabilization, 1: artificial diffusion, 2: SUPG (not working), 4: flux corrected transport"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "maxiter", "maximum number of non linear iterations"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "elementinterp", "interpolation scheme for finite elements [''P1'', ''P2'']"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "stress_threshold", "stress threshold for damage initiation (Pa)"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "stress_ubound", "stress upper bound for damage healing (Pa)"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "kappa", "ductility parameter for stress softening and damage [ > 1]"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "c1", "damage parameter 1 "))
        s += '{}\n'.format(class_utils.fielddisplay(self, "c2", "damage parameter 2 "))
        s += '{}\n'.format(class_utils.fielddisplay(self, "c3", "damage parameter 3 "))
        s += '{}\n'.format(class_utils.fielddisplay(self, "c4", "damage parameter 4 "))
        s += '{}\n'.format(class_utils.fielddisplay(self, "healing", "damage healing parameter"))
        s += '{}\n'.format(class_utils.fielddisplay(self, "equiv_stress", "0: von Mises, 1: max principal"))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - damage Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude damage fields to 3D
        """
        self.D = mesh._project_3d(md, vector = self.D, type = 'node')
        self.spcdamage = mesh._project_3d(md, vector = self.spcdamage, type = 'node')
            
        return self

    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [damage] parameters.

        Parameters
        ----------
        md : :class:`pyissm.model.Model`
            The model object to check.
        solution : :class:`pyissm.model.solution`
            The solution object to check.
        analyses : list of :class:`str`
            List of analyses to check consistency for.

        Returns 
        -------
        md : :class:`pyissm.model.Model`
            The model object with any consistency errors noted.
        """

        class_utils.check_field(md, fieldname = "damage.isdamage", scalar = True, values = [0, 1])

        if self.isdamage:
            class_utils.check_field(md, fieldname = "damage.D", ge = 0, le = self.max_damage, size = (md.mesh.numberofvertices, ))
            class_utils.check_field(md, fieldname = "damage.max_damage", ge = 0, lt = 1)
            class_utils.check_field(md, fieldname = "damage.law", scalar = True, values = [0, 1, 2, 3])
            class_utils.check_field(md, fieldname = "damage.spcdamage", allow_inf = True, timeseries = True)
            class_utils.check_field(md, fieldname = "damage.stabilization", scalar = True, values = [0, 1, 2, 4])
            class_utils.check_field(md, fieldname = "damage.maxiter", ge = 0)
            class_utils.check_field(md, fieldname = "damage.elementinterp", values = ["P1", "P2"])
            class_utils.check_field(md, fieldname = "damage.stress_threshold", ge = 0)
            class_utils.check_field(md, fieldname = "damage.stress_ubound", ge = 0)
            class_utils.check_field(md, fieldname = "damage.kappa", gt = 1)
            class_utils.check_field(md, fieldname = "damage.healing", ge = 0)
            class_utils.check_field(md, fieldname = "damage.c1", ge = 0)
            class_utils.check_field(md, fieldname = "damage.c2", ge = 0)
            class_utils.check_field(md, fieldname = "damage.c3", ge = 0)
            class_utils.check_field(md, fieldname = "damage.c4", ge = 0)
            class_utils.check_field(md, fieldname = "damage.equiv_stress", scalar = True, values = [0, 1])
            class_utils.check_field(md, fieldname = "damage.requested_outputs", string_list = True)

        elif self.law != 0:
            if solution == "DamageEvolutionSolution":
                raise RuntimeError("pyissm.model.classes.damage.check_consistency: Invalid evolution law (md.damage.law) for a damage solution")

        return md
    
    # Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self,
                        md = None,
                        return_default_outputs = False):
        """
        Process requested outputs for [damage] parameters, expanding 'default' to appropriate outputs.

        Parameters
        ----------
        md : :class:`pyissm.model.Model`, optional
            Model object containing mesh information.
        return_default_outputs : :class:`bool`, default=False
            Whether to also return the list of default outputs.
            
        Returns
        -------
        outputs : :class:`list`
            List of output strings with 'default' expanded to actual output names.
        default_outputs : :class:`list`, optional
            Returned only if `return_default_outputs=True`.
        """

        outputs = []

        ## Set default_outputs
        if md.mesh.domain_type() == '2Dhorizontal':
            default_outputs = ['DamageDbar']
        else:
            default_outputs = ['DamageD']

        ## Loop through all requested outputs
        for item in self.requested_outputs:
            
            ## Process default outputs
            if item == 'default':                
                outputs.extend(default_outputs)
            else:
                outputs.extend(default_outputs)

            ## Append other requested outputs (not defaults)
        else:
            outputs.append(item)

        if return_default_outputs:
            return outputs, default_outputs
        return outputs

    # Marshall method for saving the damage parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [damage] parameters to a binary file.

        Parameters
        ----------
        fid : :class:`file object`
            The file object to write the binary data to.
        prefix : :class:`str`
            Prefix string used for data identification in the binary file.
        md : :class:`pyissm.model.Model`, optional
            ISSM model object needed in some cases.
            
        Returns
        -------
        None
        """
        
        ## Write Boolean fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'isdamage', format = 'Boolean')

        ## If damage is enabled, write additional fields
        if self.isdamage:
            ## Write Integer fields
            execute.WriteData(fid, prefix, obj = self, fieldname = 'law', format = 'Integer')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'stabilization', format = 'Integer')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'maxiter', format = 'Integer')
            execute.WriteData(fid, prefix, obj = self, fieldname = 'equiv_stress', format = 'Integer')

            ## Write DoubleMat fields
            execute.WriteData(fid, prefix, obj = self, fieldname = 'D', format = 'DoubleMat', mattype = 1)
            execute.WriteData(fid, prefix, obj = self, fieldname = 'spcdamage', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)

            ## Write Double fields
            fieldnames = ['max_damage', 'stress_threshold', 'stress_ubound', 'kappa', 'c1', 'c2', 'c3', 'c4', 'healing']
            for fieldname in fieldnames:
                execute.WriteData(fid, prefix, obj = self, fieldname = fieldname, format = 'Double')

            ## Write other fields
            execute.WriteData(fid, prefix, name = 'md.damage.elementinterp', data = self.elementinterp, format = 'String')
            execute.WriteData(fid, prefix, name = 'md.damage.requested_outputs', data = self.process_outputs(md), format = 'StringArray')