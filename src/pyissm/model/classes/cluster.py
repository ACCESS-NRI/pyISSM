import os
import yaml
import subprocess
import numpy as np
import warnings
from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm import model, tools

## ------------------------------------------------------
## cluster.generic
## ------------------------------------------------------
@class_registry.register_class
class generic(class_registry.manage_state):
    
    """
    Generic cluster class for ISSM.

    This class provides a generic interface for managing cluster configurations
    and job execution in the Ice Sheet System Model (ISSM) framework. It handles cluster
    parameters, queue script generation, job submission, and result retrieval.
    
    Parameters
    ----------
    config_file : str, optional
        Path to YAML configuration file containing cluster parameters.
        If provided, will override default parameters with values from the file.
    other : object, optional
        Another cluster object to inherit matching fields from.
    
    Attributes
    ----------
    name : str
        Name of the cluster (defaults to hostname).
    login : str
        Login username for the cluster (defaults to current username).
    np : int
        Number of processors to use (default: 1).
    port : int
        Port number for connections (default: 0).
    interactive : int
        Interactive mode flag (default: 1).
    codepath : str
        Path to the ISSM executables directory (default: $ISSM_DIR/bin).
    executionpath : str
        Path to the execution directory on the cluster (default: $ISSM_DIR/execution).
    valgrind : str
        Path to valgrind executable for memory debugging (default: $ISSM_DIR/externalpackages/valgrind/bin/valgrind).
    valgrindlib : str
        Path to valgrind MPI debug library (default: $ISSM_DIR/externalpackages/valgrind/install/lib/libmpidebug.so).
    valgrindsup : list of str
        List of valgrind suppression files (default: $ISSM_DIR/externalpackages/valgrind/issm.supp).
    verbose : int
        Verbose output flag (default: 1).
    shell : str
        Shell to use for command execution (default: '/bin/sh').

    Methods
    -------
    build_queue_script(dir_name, model_name, solution, io_gather, is_valgrind,
                       is_gprof, is_dakota, is_ocean_coupling, executable='issm.exe')
        Generate queue script for model execution.
    build_kriging_queue_script(model_name, solution, io_gather, is_valgrind,
                               is_gprof, executable='kriging.exe')
        Generate queue script for kriging execution.
    upload_queue_job(model_name, dir_name, file_list)
        Upload job files to the cluster.
    launch_queue_job(model_name, dir_name, restart=None, batch=False)
        Launch job execution on the cluster.
    download(dir_name, file_list)
        Download job results from the cluster.
    
    Notes
    -----
    This class inherits from class_registry.manage_state to provide state
    management capabilities. Configuration parameters can be overridden
    via YAML configuration files or by inheriting from other cluster objects.

    Examples
    --------
    >>> cluster = generic()
    >>> cluster.np = 4
    >>> cluster.name = 'my_cluster'
    """

    # Initialise with default parameters
    def __init__(self, config_file = None, other = None):
        self.name = tools.config.get_hostname()
        self.login = tools.config.get_username()
        self.np = 1
        self.port = 0
        self.interactive = 1

        # Set paths based on ISSM_DIR
        if tools.config.get_issm_dir() is None:
            warnings.warn('pyissm.model.classes.cluster.generic: ISSM_DIR is not set. Returning empty cluster paths.\n')
            self.codepath = ''
            self.executionpath = ''
            self.valgrind = ''
            self.valgrindlib = ''
            self.valgrindsup = []
        else:
            self.codepath = tools.config.get_issm_dir() + '/bin'
            self.executionpath = tools.config.get_issm_dir() + '/execution'
            self.valgrind = tools.config.get_issm_dir() + '/externalpackages/valgrind/bin/valgrind'
            self.valgrindlib = tools.config.get_issm_dir() + '/externalpackages/valgrind/install/lib/libmpidebug.so'
            self.valgrindsup = [tools.config.get_issm_dir() + '/externalpackages/valgrind/issm.supp']  # add any .supp in list form as needed
        
        self.verbose = 1
        self.shell = '/bin/sh'

        # Inherit matching fields from provided class
        super().__init__(other)

        # Override default parameters with config file values
        if config_file is not None:
            ## Load yaml file
            with open(config_file, 'r') as fid:
                cfg = yaml.safe_load(fid)

            ## Set attributes (Ignore keys that do not match any attribute)
            for key, value in cfg.items():
                if hasattr(self, key):
                    setattr(self, key, value)

    # Define repr
    def __repr__(self):
        s = '   Cluster parameters:\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'name', 'name of the cluster'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'login', 'login name for the cluster'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'np', 'number of processors'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'port', 'port number'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'interactive', 'interactive mode'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'codepath', 'path to the code'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'executionpath', 'path to the execution directory'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'valgrind', 'path to valgrind'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'valgrindlib', 'path to valgrind library'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'valgrindsup', 'valgrind suppression files'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'verbose', 'verbose mode'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'shell', 'shell to use'))

        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - cluster.generic Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        if self.np < 1 or np.isnan(self.np):
            md.check_message('cluster.generic.check_consistency: np must be >= 1')

        pwd = os.getcwd()
        if self.executionpath == pwd:
            md.check_message('cluster.generic.check_consistency: md.cluster.executionpath must be different from the current working directory.')

        return md

    def build_queue_script(self,
                           dir_name,
                           model_name,
                           solution,
                           io_gather,
                           is_valgrind,
                           is_gprof,
                           is_dakota,
                           is_ocean_coupling,
                           executable = 'issm.exe'):
        """
        Build a queue script for executing ISSM models on the cluster.
        
        This method generates platform-specific execution scripts (bash for Linux/Mac, 
        batch for Windows) that handle model execution with various configurations 
        including MPI, debugging tools, and specialized executables.
        Parameters
        ----------
        dir_name : str
            Directory name where the model files are located.
        model_name : str
            Name of the model to execute.
        solution : str
            Solution type or configuration parameter.
        io_gather : bool
            Flag indicating whether to gather I/O operations. If False, output 
            files will be concatenated.
        is_valgrind : bool
            Flag to enable Valgrind memory debugging tool execution.
        is_gprof : bool
            Flag to enable gprof profiling tool execution.
        is_dakota : bool
            Flag to use DAKOTA optimization executable.
        is_ocean_coupling : bool
            Flag to use ocean coupling executable.
        executable : str, optional
            Name of the executable file to run. Default is 'issm.exe'.

        Raises
        ------
        IOError
            If Python wrappers are not installed or if DAKOTA support is 
            requested but not available in the ISSM build.

        Notes
        -----
        - On Linux/Mac systems, creates a '.queue' bash script
        - On Windows systems, creates a '.bat' batch script
        - Automatically handles MPI execution when available
        - In interactive mode, creates empty error and output log files
        - Supports various debugging and profiling tools integration
        - Handles different executable types based on coupling requirements
        """
        
        # Require wrappers when executing a model
        if not tools.wrappers.check_wrappers_installed():
            raise IOError('cluster.generic.build_queue_script: Python wrappers not installed. Unable to build queue script.')

        # DAKOTA EXECUTABLE
        if is_dakota:
            ## Check that ISSM has DAKOTA support and get version #
            if tools.wrappers.IssmConfig('_HAVE_DAKOTA_')[0]:
                version = float(tools.wrappers.IssmConfig('_DAKOTA_VERSION_')[0])
                if version >= 6:
                    executable = 'issm_dakota.exe'
            else:
                ## If no dakota support, raise error
                raise IOError('cluster.generic.build_queue_script error: ISSM not built with DAKOTA support')
        
        # OCEAN COUPLING EXECUTABLE
        if is_ocean_coupling:
            executable = 'issm_ocean_coupling.exe'

        # BUILD SCRIPT
        ## Linux/Mac
        if not tools.config.is_pc():
            fid = open(model_name + '.queue', 'w')
            fid.write('#!/bin/bash\n')

            if not is_valgrind:
                if self.interactive:
                    if tools.wrappers.IssmConfig('_HAVE_MPI_')[0]:
                        fid.write('mpiexec -np {} {}/{} {} {}/{} {}'.format(self.np, self.codepath, executable, solution, self.executionpath, dir_name, model_name))
                    else:
                        fid.write('{}/{} {} {}/{} {}'.format(self.codepath, executable, solution, self.executionpath, dir_name, model_name))
                else:
                    if tools.wrappers.IssmConfig('_HAVE_MPI_')[0]:
                        fid.write('mpiexec -np {} {}/{} {} {}/{} {} 2> {}.errlog > {}.outlog'.format(self.np, self.codepath, executable, solution, self.executionpath, dir_name, model_name, model_name, model_name))
                    else:
                        fid.write('{}/{} {} {}/{} {} 2> {}.errlog > {}.outlog '.format(self.codepath, executable, solution, self.executionpath, dir_name, model_name, model_name, model_name))
            elif is_gprof:
                fid.write('\n gprof {}/{} gmon.out > {}.performance'.format(self.codepath, executable, model_name))
            else:
                supstring = ''
                for supfile in self.valgrindsup:
                    supstring += ' --suppressions=' + supfile
                
                if tools.wrappers.IssmConfig('_HAVE_MPI_')[0]:
                    fid.write('mpiexec -np {} {} --leak-check=full {} {}/{} {} {}/{} {} 2> {}.errlog > {}.outlog '.format(self.np, self.valgrind, supstring, self.codepath, executable, solution, self.executionpath, dir_name, model_name, model_name, model_name))
                else:
                    fid.write('{} --leak-check=full {} {}/{} {} {}/{} {} 2> {}.errlog > {}.outlog '.format(self.valgrind, supstring, self.codepath, executable, solution, self.executionpath, dir_name, model_name, model_name, model_name))

            # Concatenate output files if needed
            if not io_gather:
                fid.write('\ncat {}.outbin .*>{}.outbin'.format(model_name, model_name))
            fid.close()
        
        ## Windows
        else:
            fid = open(model_name + '.bat', 'w')
            fid.write('@echo off\n')
            if self.interactive:
                fid.write('"{}/{}" {} "{}/{}" {} '.format(self.codepath, executable, solution, self.executionpath, dir_name, model_name))
            else:
                fid.write('"{}/{}" {} "{}/{}" {} 2>{}.errlog>{}.outlog'.format(self.codepath, executable, solution, self.executionpath, dir_name, model_name, model_name, model_name))
            fid.close()

        ## In interactive mode, create a run file, and errlog and outlog file
        if self.interactive:
            fid = open(model_name + '.errlog', 'w')
            fid.close()
            fid = open(model_name + '.outlog', 'w')
            fid.close()

    def build_kriging_queue_script(self,
                                   model_name,
                                   solution,
                                   io_gather,
                                   is_valgrind,
                                   is_gprof,
                                   executable = 'kriging.exe'):
        
        """
        Build a queue script for executing kriging models on the cluster.
        
        This method generates platform-specific execution scripts (bash for Linux/Mac, 
        batch for Windows) that handle kriging model execution with various configurations 
        including MPI, debugging tools, and profiling.
        
        Parameters
        ----------
        model_name : str
            Name of the kriging model to execute.
        solution : str
            Solution type or configuration parameter.
        io_gather : bool
            Flag indicating whether to gather I/O operations. If False, output 
            files will be concatenated.
        is_valgrind : bool
            Flag to enable Valgrind memory debugging tool execution.
        is_gprof : bool
            Flag to enable gprof profiling tool execution.
        executable : str, optional
            Name of the executable file to run. Default is 'kriging.exe'.
        
        Raises
        ------
        IOError
            If Python wrappers are not installed.
        
        Notes
        -----
        - On Linux/Mac systems, creates a '.queue' bash script
        - On Windows systems, creates a '.bat' batch script
        - Automatically handles MPI execution for kriging operations
        - In interactive mode, creates empty error and output log files
        - Supports memory debugging with Valgrind and profiling with gprof
        - Specifically designed for kriging executable execution
        """
        
        # Require wrappers when executing a model
        if not tools.wrappers.check_wrappers_installed():
            raise IOError('cluster.generic.build_kriging_queue_script: Python wrappers not installed. Unable to build queue script.')
        
        # BUILD SCRIPT
        ## Linux/Mac
        if not tools.config.is_pc():
            fid = open(model_name + '.queue', 'w')
            fid.write('#!/bin/bash\n')
        
            if not is_valgrind:
                if self.interactive:
                    fid.write('mpiexec -np {} {}/{} {}/{} {} '.format(self.np, self.codepath, executable, self.executionpath, model_name, model_name))
                else:
                    fid.write('mpiexec -np {} {}/{} {}/{} {} 2>{}.errlog>{}.outlog '.format(self.np, self.codepath, executable, self.executionpath, model_name, model_name, model_name, model_name))
            elif is_gprof:
                fid.write('\n gprof {}/{} gmon.out>{}.performance'.format(self.codepath, executable, model_name))
            else:
                #Add --    gen - suppressions = all to get suppression lines
                #fid.write('LD_PRELOAD={} \\\n'.format(self.valgrindlib))
                fid.write('mpiexec -np {} {} --leak -check=full --suppressions={} {}/{}} {}/{} {} 2 > {}.errlog > {}.outlog ' .format(self.np, self.valgrind, self.valgrindsup, self.codepath, executable, self.executionpath, model_name, model_name, model_name, model_name))
            if not io_gather:    #concatenate the output files:
                fid.write('\ncat {}.outbin. *>{}.outbin'.format(model_name, model_name))
            fid.close()
        
        ## Windows
        else:
            fid = open(model_name + '.bat', 'w')
            fid.write('@echo off\n')
            if self.interactive:
                fid.write('"{}/issm.exe" {} "{}/{}" {} '.format(self.codepath, solution, self.executionpath, model_name, model_name))
            else:
                fid.write('"{}/issm.exe" {} "{}/{}" {} 2>{}.errlog>{}.outlog'.format
                          (self.codepath, solution, self.executionpath, model_name, model_name, model_name, model_name))
            fid.close()

        ## In interactive mode, create a run file, and errlog and outlog file
        if self.interactive:
            fid = open(model_name + '.errlog', 'w')
            fid.close()
            fid = open(model_name + '.outlog', 'w')
            fid.close()

    def upload_queue_job(self, model_name, dir_name, file_list):
        """
        Upload job files to the cluster queue system.

        Compresses the specified files into a tar.gz archive and transfers it to the
        cluster using SCP. If running in interactive mode, also includes error and
        output log files in the archive.
        
        Parameters
        ----------
        model_name : str
            Name of the model, used for naming log files in interactive mode.
        dir_name : str
            Name of the directory/archive to be created (without extension).
        file_list : list of str
            List of file paths to be included in the compressed archive.
        
        Notes
        -----
        The function creates a tar.gz archive with the name `{dir_name}.tar.gz`
        containing all files in `file_list`. In interactive mode, it also includes
        `{model_name}.errlog` and `{model_name}.outlog` files.
        The compressed archive is then transferred to the cluster using the
        cluster's configured connection parameters (name, execution path, login,
        and port).
        
        See Also
        --------
        model.io.issm_scp_out : Function used for transferring files to cluster
        """

        ## Compress files into one zip
        compress_string = 'tar -zcf {}.tar.gz '.format(dir_name)

        for file in file_list:
            compress_string += ' {} '.format(file)
        if self.interactive:
            compress_string += ' {}.errlog {}.outlog'.format(model_name, model_name)
        subprocess.call(compress_string, shell = True)

        ## Transfer to cluster
        print(f'Transferring {dir_name}.tar.gz to cluster {self.name}...')
        model.io.issm_scp_out(self.name, self.executionpath, self.login, self.port, [dir_name + '.tar.gz'])
        
    def launch_queue_job(self,
                         model_name,
                         dir_name,
                         restart = None,
                         batch = False):
        """
        Launch a job on the cluster queue system.
        This method builds and executes the appropriate launch command for submitting
        a job to the cluster's queue system. It handles both fresh job submissions
        and job restarts, with optional batch processing mode.
        Parameters
        ----------
        model_name : str
            Name of the model to be executed on the cluster.
        dir_name : str
            Name of the directory where the job will be executed.
        restart : bool or None, optional
            If not None, indicates this is a restart of an existing job.
            When restarting, the method assumes the job directory already exists
            and only executes the queue script. Default is None.
        batch : bool, optional
            Flag indicating whether to run in batch mode. When True, only
            extracts the tar.gz file without executing the queue script.
            When False (default), extracts and immediately executes the job.
            Only relevant when restart is None.
        Notes
        -----
        The method performs different operations based on the restart parameter:
        - If restart is not None: Changes to the execution directory and runs
            the existing queue script.
        - If restart is None: Removes any existing directory, creates a new one,
            moves and extracts the tar.gz file, and optionally runs the queue script
            depending on the batch parameter.
        The job is launched via SSH connection to the cluster using the cluster's
        name, login credentials, and port configuration.
        Examples
        --------
        Launch a new job:
        >>> cluster.launch_queue_job('simulation_01', 'run_dir')
        Restart an existing job:
        >>> cluster.launch_queue_job('simulation_01', 'run_dir', restart=True)
        Launch in batch mode (extract only, no execution):
        >>> cluster.launch_queue_job('simulation_01', 'run_dir', batch=True)
        """

        # Build launch command        
        ## If not restarting, move to execution path and execute model
        if restart is not None:
            launch_command = 'cd {} && cd {} chmod 755 {}.queue && ./{}.queue'.format(self.executionpath, dir_name, model_name, model_name)

        ## If not restarting, remove existing directory, recreate it, move tar.gz file and uncompress
        else:
            if batch:
                launch_command = 'cd {} && rm -rf ./{} && mkdir {} && cd {} && mv ../{}.tar.gz ./&& tar -zxf {}.tar.gz'.format(self.executionpath, dir_name, dir_name, dir_name, dir_name, dir_name)
            else:
                launch_command = 'cd {} && rm -rf ./{} && mkdir {} && cd {} && mv ../{}.tar.gz ./&& tar -zxf {}.tar.gz  && chmod 755 {}.queue && ./{}.queue'.format(self.executionpath, dir_name, dir_name, dir_name, dir_name, dir_name, model_name, model_name)
                   
        ## Launch job on cluster
        print(f'Launching job {model_name} on cluster {self.name}...')
        model.io.issm_ssh(self.name, self.login, self.port, launch_command)

    def download(self, dir_name, file_list):
        """
        Download files from a remote cluster to the local machine.
        This method retrieves specified files from a remote cluster directory
        to the current local directory. On Windows systems, this operation
        is skipped as it's not supported.

        Parameters
        ----------
        dir_name : str
            The name of the directory on the remote cluster containing the files
            to download.
        file_list : list of str
            A list of filenames to download from the remote cluster directory.

        Returns
        -------
        None

        Notes
        -----
        - This method does nothing on Windows platforms and returns immediately.
        - Files are copied from the cluster's execution path combined with the
          specified directory name.
        - The actual file transfer is handled by the `model.io.issm_scp_in` function.
        """

        if tools.config.is_pc():
            ## Do nothing on Windows
            return
        
        # Copy files from cluster to current directory
        print(f'Retrieving results from cluster {self.name}...')
        directory = f'{self.executionpath}/{dir_name}/'
        model.io.issm_scp_in(self.name, self.login, self.port, directory, file_list)

## ------------------------------------------------------
## cluster.gadi
## ------------------------------------------------------
@class_registry.register_class
class gadi(class_registry.manage_state):
    """
    Gadi HPC cluster interface for ISSM job submission and management.

    This class represents the Gadi HPC cluster at the National Computational 
    Infrastructure (NCI) and provides methods for configuring cluster parameters, 
    building PBS queue scripts, and managing job submission and file transfers.

    The Gadi cluster uses PBS Pro for job scheduling and supports parallel 
    execution via MPI. Configuration can be provided via YAML config files or 
    programmatically through object attributes.

    Parameters
    ----------
    config_file : str, optional
        Path to YAML configuration file containing cluster parameters.
        If provided, will override default parameters with values from the file.
    other : object, optional
        Another cluster object to inherit matching fields from.

    Attributes
    ----------
    name : str
        Hostname of the cluster. Defaults to 'gadi.nci.org.au' if not on Gadi.
    login : str
        Login username for the cluster. Must be provided for cluster access.
    np : int
        Number of processors to use for job execution. Default is 16.
    memory : int
        Memory per node in GB. Default is 40.
    port : int
        SSH port number for cluster connection. Default is 0.
    queue : str
        PBS queue name. Options include 'normal', 'express', 'hugemem'. 
        Default is 'normal'.
    time : int
        Walltime limit for job execution in minutes. Default is 60.
    codepath : str
        Path to the ISSM executable directory (e.g., $ISSM_DIR/bin). 
        Must be provided.
    executionpath : str
        Path to the execution/working directory on the cluster. Must be provided.
    project : str
        NCI project code for job submission. Must be provided.
    storage : str
        Storage paths to access (e.g., 'gdata/XXX+scratch/XXX'). 
        Must be provided.
    moduleload : list of str
        List of module load commands needed for PBS job execution.
    moduleuse : list of str
        List of module use commands to specify module paths.

    Notes
    -----
    All required attributes (login, codepath, executionpath, project, storage) 
    must be set either via configuration file or programmatically before building 
    and launching jobs. The moduleload and moduleuse lists must have equal length.

    Queue specifications:
    
    - normal: 48 hours on up to 3072 cores
    - express: 2 hours on up to 960 cores
    - hugemem: 48 hours on up to 3072 cores

    Examples
    --------
    >>> cluster = gadi(config_file='gadi_config.yaml')
    >>> cluster.np = 32
    >>> cluster.queue = 'express'
    """

    # Initialise with default parameters
    def __init__(self, config_file = None, other = None):

        # Adjust hostname for use on/off gadi
        host_name = tools.config.get_hostname()
        if 'gadi' not in host_name:
            host_name = 'gadi.nci.org.au'

        # Define default parameters
        self.name = host_name
        self.login = ''
        self.np = 16
        self.memory = 40
        self.port = 0
        self.queue = 'normal'
        self.time = 60 
        self.codepath = ''
        self.executionpath = ''
        self.project = ''
        self.storage = ''
        self.moduleload = []
        self.moduleuse = []

        # Inherit matching fields from provided class
        super().__init__(other)

        # Override default parameters with config file values
        if config_file is not None:
            ## Load yaml file
            with open(config_file, 'r') as fid:
                cfg = yaml.safe_load(fid)

            ## Set attributes (Ignore keys that do not match any attribute)
            for key, value in cfg.items():
                if hasattr(self, key):
                    setattr(self, key, value)
        
    # Define repr
    def __repr__(self):
        s = '   Cluster parameters:\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'name', 'name of the cluster'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'login', 'login name for the cluster'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'np', 'number of processors'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'memory', 'memory per node (in GB)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'port', 'port number'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'queue', 'queue name'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'time', 'walltime (in minutes)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'codepath', 'path to the ISSM executable (e.g. $ISSM_DIR/bin)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'executionpath', 'path to the execution directory'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'project', 'NCI project name'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'storage', 'storage path (e.g. gdata/XXX+scratch/XXX)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'moduleload', 'List of module load commands needed for PBS job'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'moduleuse', 'List of module use commands needed for PBS job'))

        return s
    
    # Define class string
    def __str__(self):
        s = 'ISSM - cluster.gadi Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):

        queue_dict = {
            'normal': [48*60, 3072], # 48h on 3072 cores
            'express': [2*60, 960], # 2h on 960 cores
            'hugemem': [48*60, 3072], # 48h on 3072 cores
        }

        class_utils.cluster_queue_requirements(queue_dict, self.queue, self.np, self.time)

        if not self.login:
            md.check_message('pyissm.model.classes.cluster.gadi: login name must be provided for gadi cluster')

        if not self.codepath:
            md.check_message('pyissm.model.classes.cluster.gadi: codepath must be provided for gadi cluster')
        
        if not self.executionpath:
            md.check_message('pyissm.model.classes.cluster.gadi: executionpath must be provided for gadi cluster')

        if not self.project:
            md.check_message('pyissm.model.classes.cluster.gadi: project must be provided for gadi cluster')

        if not self.storage:
            md.check_message('pyissm.model.classes.cluster.gadi: storage must be provided for gadi cluster')

        if len(self.moduleload) != len(self.moduleuse):
            md.check_message('pyissm.model.classes.cluster.gadi: moduleload and moduleuse must have the same length')

        pwd = os.getcwd()
        if self.executionpath == pwd:
            md.check_message('cluster.generic.check_consistency: md.cluster.executionpath must be different from the current working directory.')


        return md

    # Build queue script
    def build_queue_script(self,
                           dir_name,
                           model_name,
                           solution,
                           io_gather,
                           is_valgrind,
                           is_gprof,
                           is_dakota,
                           is_ocean_coupling,
                           executable = 'issm.exe'):
        """
        Generate a PBS queue submission script for running ISSM models
        on the Gadi cluster. The script includes resource specifications, module
        loading, and execution commands.

        Parameters
        ----------
        dir_name : str
            Directory name where the model execution files are stored.
        model_name : str
            Name of the model, used for output file naming.
        solution : str
            Solution type or identifier to pass to the executable.
        io_gather : bool
            If True, output files are pre-gathered. If False, output binary files
            are concatenated after execution.
        is_valgrind : bool
            If True, raises NotImplementedError as Valgrind is not supported.
        is_gprof : bool
            If True, raises NotImplementedError as gprof is not supported.
        is_dakota : bool
            If True, raises NotImplementedError as DAKOTA is not supported.
        is_ocean_coupling : bool
            If True, raises NotImplementedError as ocean coupling is not supported.
        executable : str, optional
            Name of the executable to run. Default is 'issm.exe'.

        Raises
        ------
        IOError
            If Python wrappers are not installed.
        NotImplementedError
            If any of the unsupported features (DAKOTA, ocean coupling, Valgrind, gprof)
            are requested.

        Returns
        -------
        None
            Writes a queue script file named '{model_name}.queue' to the current directory.
        """
        
        # Require wrappers when executing a model
        if not tools.wrappers.check_wrappers_installed():
            raise IOError('pyissm.model.classes.cluster.gadi.build_queue_script: Python wrappers not installed. Unable to build queue script.')

        # Raise error for unimplemented features
        if is_dakota:
            raise NotImplementedError('pyissm.model.classes.cluster.gadi.build_queue_script: DAKOTA support not implemented for gadi cluster yet.')
        
        if is_ocean_coupling:
            raise NotImplementedError('pyissm.model.classes.cluster.gadi.build_queue_script: Ocean coupling support not implemented for gadi cluster yet.')
        
        if is_valgrind:
            raise NotImplementedError('pyissm.model.classes.cluster.gadi.build_queue_script: Valgrind support not implemented for gadi cluster yet.')
        
        if is_gprof:
            raise NotImplementedError('pyissm.model.classes.cluster.gadi.build_queue_script: gprof support not implemented for gadi cluster yet.')
        
        # Write queue script
        fid = open(model_name + '.queue', 'w')
        fid.write('#!/bin/bash\n')
        fid.write(f'#PBS -P {self.project}\n')
        fid.write(f'#PBS -q {self.queue}\n')
        fid.write(f'#PBS -l ncpus={self.np}\n')
        fid.write(f'#PBS -l mem={self.memory}GB\n')
        fid.write(f'#PBS -l walltime={self.time*60}\n') # Walltime is in seconds
        fid.write('#PBS -l wd\n')  
        fid.write('#PBS -j oe\n')
        fid.write(f'#PBS -l storage={self.storage}\n')
        fid.write(f'#PBS -o {self.executionpath}/{dir_name}/{model_name}.outlog \n')
        fid.write(f'#PBS -e {self.executionpath}/{dir_name}/{model_name}.errlog \n\n')   

        # Load modules as needed (print alternating module use and module load lines)
        fid.write('module purge\n')
        for x, y in zip(self.moduleuse, self.moduleload):
            fid.write(f"module use {x}\n")
            fid.write(f"module load {y}\n")

        # Add execution command
        fid.write(f'mpiexec -n {self.np} {self.codepath}/{executable} {solution} {self.executionpath}/{dir_name} {model_name}\n')

        # Concatenate output files if needed
        if not io_gather:
            fid.write(f'\ncat {model_name}.outbin.* > {model_name}.outbin\n')
        # Close file
        fid.close()

    # Upload job to cluster
    def upload_queue_job(self, model_name, dir_name, file_list):
        """
        Upload job files to the cluster queue system.

        This method compresses the specified files into a tar.gz archive and transfers it to the
        cluster using SCP. If running in interactive mode, also includes error and
        output log files in the archive.

        Parameters
        ----------
        model_name : str
            Name of the model, used for naming log files in interactive mode (not used here).
        dir_name : str
            Name of the directory/archive to be created (without extension).
        file_list : list of str
            List of file paths to be included in the compressed archive.

        Notes
        -----
        The function creates a tar.gz archive with the name `{dir_name}.tar.gz`
        containing all files in `file_list`. The compressed archive is then transferred to the cluster using the
        cluster's configured connection parameters (name, execution path, login,
        and port).
        """

        # Compress inputs into a tarball
        compressstring = f'tar -zcf {dir_name}.tar.gz'
        for f in file_list:
            compressstring += f' {f}'

        subprocess.call(compressstring, shell = True)
        print(f'Transferring {dir_name}.tar.gz to cluster {self.name}...')
        directory = self.executionpath
        model.io.issm_scp_out(self.name, directory, self.login, self.port, [dir_name + '.tar.gz'])

    # Launch job on cluster
    def launch_queue_job(self,
                         model_name,
                         dir_name,
                         restart = None,
                         batch = False):
        """
        Launch a job on the Gadi cluster queue system.
        
        This method submits a job to the Gadi PBS queue system. It handles both
        fresh job submissions and job restarts, with optional batch processing mode.

        Parameters
        ----------
        model_name : str
            Name of the model to be executed on the cluster.
        dir_name : str
            Name of the directory where the job will be executed.
        restart : bool or None, optional
            If not None, indicates this is a restart of an existing job.
            When restarting, the method assumes the job directory already exists
            and only submits the queue script via qsub. Default is None.
        batch : bool, optional
            Flag indicating whether to run in batch mode. Currently unused for
            Gadi cluster but maintained for interface compatibility.
            Default is False.

        Returns
        -------
        None

        Notes
        -----
        The method performs different operations based on the restart parameter:
        
        - If restart is not None: Changes to the existing execution directory and
          submits the queue script using qsub.
        - If restart is None: Removes any existing directory, creates a new one,
          moves and extracts the tar.gz file, then submits the queue script using qsub.
        
        The job is launched via SSH connection to the cluster using the cluster's
        name, login credentials, and port configuration.

        Examples
        --------
        Launch a new job:
        
        >>> cluster.launch_queue_job('simulation_01', 'run_dir')
        
        Restart an existing job:
        
        >>> cluster.launch_queue_job('simulation_01', 'run_dir', restart=True)
        """
        
        if restart is not None:
            # Just qsub in existing directory
            launch_command = (
                f'cd {self.executionpath}/{dir_name} && qsub {model_name}.queue')
        else:
            # Create/clean directory, extract tar, then qsub
            launch_command = (
                f'cd {self.executionpath} && '
                f'rm -rf ./{dir_name} && mkdir {dir_name} && cd {dir_name} && '
                f'mv ../{dir_name}.tar.gz ./ && '
                f'tar -zxf {dir_name}.tar.gz && '
                f'qsub {model_name}.queue'
                )
            
        print(f'Launching job {model_name} on cluster {self.name}...')
        model.io.issm_ssh(self.name, self.login, self.port, launch_command)

    # Download results from cluster
    def download(self, dir_name, file_list):
        """
        Download files from a remote cluster to the local machine.

        Parameters
        ----------
        dir_name : str
            The name of the directory on the remote cluster containing the files
            to download.
        file_list : list of str
            A list of filenames to download from the remote cluster directory.

        Returns
        -------
        None
        """

        # Copy files from cluster to current directory
        print(f'Retrieving results from cluster {self.name}...')
        directory = f'{self.executionpath}/{dir_name}/'
        model.io.issm_scp_in(self.name, self.login, self.port, directory, file_list)