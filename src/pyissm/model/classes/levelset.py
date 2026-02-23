import numpy as np
from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm.model import execute, mesh

@class_registry.register_class
class levelset(class_registry.manage_state):
    """
    Level-set method parameters class for ISSM.

    This class encapsulates parameters for configuring level-set method simulations in the ISSM (Ice Sheet System Model) framework.
    The level-set method is used to track moving boundaries such as ice fronts and calving fronts, 
    allowing for dynamic changes in ice sheet geometry during simulations.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    stabilization : int, default=1
        Stabilization method: 0=No Stabilization, 1=Artificial Diffusivity (most stable, least accurate), 2=Streamline Upwinding, 5=SUPG (most accurate, may be unstable).
    spclevelset : ndarray, default=nan
        Levelset constraints (NaN means no constraint).
    reinit_frequency : int, default=10
        Amount of time steps after which the levelset function is re-initialized.
    kill_icebergs : int, default=1
        Remove floating icebergs to prevent rigid body motions (1=true, 0=false).
    migration_max : float, default=1e12
        Maximum allowed migration rate [m/a].
    fe : str, default='P1'
        Finite Element type: 'P1' (default) or 'P2'.

    Methods
    -------
    __init__(self, other=None)
        Initializes the levelset parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the levelset parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.levelset = pyissm.model.classes.levelset()
    md.levelset.stabilization = 5
    md.levelset.reinit_frequency = 5
    md.levelset.kill_icebergs = 0
    md.levelset.migration_max = 1e10
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.stabilization = 1
        self.spclevelset = np.nan
        self.reinit_frequency = 10
        self.kill_icebergs = 1
        self.migration_max = 1e12
        self.fe = 'P1'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Level-set parameters:\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'stabilization', '0: No Stabilization - No stabilization techniques applied.'))
        s += '{}\n'.format('                             1: Artificial Diffusivity - Most stable, but least accurate.')
        s += '{}\n'.format('                             2: Streamline Upwinding')
        s += '{}\n'.format('                             5: SUPG - Most accurate, but may be unstable in some applications.')
        s += '{}\n'.format(class_utils.fielddisplay(self, 'spclevelset', 'Levelset constraints (NaN means no constraint)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'reinit_frequency', 'Amount of time steps after which the levelset function in re-initialized'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'kill_icebergs', 'remove floating icebergs to prevent rigid body motions (1: true, 0: false)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'migration_max', 'maximum allowed migration rate (m/a)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'fe', 'Finite Element type: \'P1\' (default), or \'P2\''))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - levelset Class'
        return s
    
    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude levelset fields to 3D
        """
        self.spclevelset = mesh.project_3d(md, vector = self.spclevelset, type = 'node')
            
        return self
    



    def reinitialize(self, md, levelset):

        import os
        import tempfile
        import numpy as np
        from pyissm.model import mesh
        from pyissm.tools.wrappers import ExpToLevelSet
        from pyissm.tools.exp import exp_write, isoline
        if levelset is None or len(levelset) == 0:
            raise RuntimeError("levelset provided is empty")

        levelset = np.asarray(levelset).reshape(-1)
        dim = md.mesh.dimension()

        if dim == 3:
            if levelset.size != md.mesh.numberofvertices2d:
                raise RuntimeError("levelset provided should be specified at the 2d vertices of the mesh")
            x = np.asarray(md.mesh.x2d).reshape(-1)
            y = np.asarray(md.mesh.y2d).reshape(-1)
        else:
            if levelset.size != md.mesh.numberofvertices:
                raise RuntimeError("levelset provided should be specified at the vertices of the mesh")
            x = np.asarray(md.mesh.x).reshape(-1)
            y = np.asarray(md.mesh.y).reshape(-1)

        # 1) extract 0-contour as python dict/list
        contours, _ = isoline(md, levelset, value=0.0)

        if not isinstance(contours, list) or len(contours) == 0:
            raise RuntimeError("isoline returned no contours at value=0")

        # 2) force .exp write so ExpToLevelSet backend can load it
        with tempfile.TemporaryDirectory() as td:
            exp_path = os.path.join(td, "levelset0.exp")
            exp_write(contours, exp_path)
            dist = ExpToLevelSet(x, y, exp_path)

        dist = np.asarray(dist).reshape(-1)

        # If ExpToLevelSet returns N+1 metadata row, drop it
        n = x.size
        dist = dist[:n]

        levelsetnew = np.abs(dist)

        # 3) restore sign
        pos = np.where(levelset < 0)[0]
        if dim == 3:
            # Expand distance to 3D if needed (depends on what ExpToLevelSet returned)
            nv2d = md.mesh.numberofvertices2d
            nl = md.mesh.numberoflayers
            # If levelsetnew is 2D length, tile to 3D length
            if levelsetnew.size == nv2d:
                levelsetnew = np.tile(levelsetnew, nl)
            for i in range(nl):
                pos3d = pos + i * nv2d
                levelsetnew[pos3d] = -levelsetnew[pos3d]
        else:
            levelsetnew[pos] = -levelsetnew[pos]

        return levelsetnew

    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        # Early return if not a transient moving front simulation
        if (solution != 'TransientSolution') or (not md.transient.ismovingfront):
            return md

        class_utils.check_field(md, fieldname = 'levelset.spclevelset', timeseries = True, allow_inf = False)
        class_utils.check_field(md, fieldname = 'levelset.stabilization', scalar = True, values = [0, 1, 2, 5, 6])
        class_utils.check_field(md, fieldname = 'levelset.kill_icebergs', scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = 'levelset.migration_max', scalar = True, gt = 0, allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = 'levelset.fe', values = ['P1', 'P2'])

        return md

    # Marshall method for saving the levelset parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [levelset] parameters to a binary file.

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
        execute.WriteData(fid, prefix, obj = self, fieldname = 'stabilization', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'spclevelset', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'reinit_frequency', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'kill_icebergs', format = 'Boolean')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'migration_max', format = 'Double', scale = 1. / md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'fe', format = 'String')
