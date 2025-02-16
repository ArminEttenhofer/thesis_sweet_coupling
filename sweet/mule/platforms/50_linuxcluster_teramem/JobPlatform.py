import platform
"""
Platform job scripts for Linux Cluster TeraMem nodes, see also

https://doku.lrz.de/large-memory-teramem-11484374.html
"""

import socket
import sys
import os

from mule.JobGeneration import *
from mule.JobPlatformResources import *
from . import JobPlatformAutodetect

# Underscore defines symbols to be private
_job_id = None



def get_platform_autodetect():
    """
    Returns
    -------
    bool
    	True if current platform matches, otherwise False
    """

    return JobPlatformAutodetect.autodetect()



def get_platform_id():
    """
    Return platform ID

    Returns
    -------
    string
    	unique ID of platform
    """

    return "linuxcluster_teramem"



def get_platform_resources():
    """
    Return information about hardware
    """

    h = JobPlatformResources()

    """
    4 x Intel Xeon Platinum 8360HL
    """

    h.num_cores_per_node = 96
    # Number of nodes per job are limited
    h.num_nodes = 1
    #h.num_nodes = 60
    h.num_cores_per_socket = 24
    h.max_wallclock_seconds = 48*60*60
    return h



def jobscript_setup(jg : JobGeneration):
    """
    Setup data to generate job script
    """

    pass



def jobscript_get_header(jg : JobGeneration):
    """
    These headers typically contain the information on e.g. Job exection, number of compute nodes, etc.

    Returns
    -------
    string
    	multiline text for scripts
    """
    _job_id = jg.runtime.getUniqueID(jg.compile, jg.unique_id_filter)

    p = jg.parallelization
    time_str = p.get_max_wallclock_seconds_hh_mm_ss()

    #
    # Next, we load the linux cluster specific configuration file which we assume to be in the HOME folder of the user
    #
    config_file = os.environ["HOME"]+"/sweet_platform_user_specific_config.py"
    config_per_user = {}

    # Special override for unit test
    if "MULE_TEST_PLATFORMS" in os.environ:
        config_per_user["project_id"] = "dummy_project_id"
        config_per_user["user_email"] = "dummy_email"
    else:
        try:
            with open(config_file, "rb") as f:
                code = compile(f.read(), config_file, "exec")
                exec(code, config_per_user)

        except:
            print("*"*80)
            print("ERROR: Failed to parse ""+config_file+""")
            print("*"*80)
            print("""

Make sure that it"s in the following format:

user_email = "someemail@address.com"
project_id = "pro123abc"        # Not required, yet - can be left to empty string

""")
            raise Exception("ERROR - stopping here")

    print("Project ID: "+config_per_user["project_id"])
    print("User Email: "+config_per_user["user_email"])

    #
    # See https://www.lrz.de/services/compute/linux-cluster/batch_parallel/example_jobs/
    #
    content = f"""#! /bin/bash
#SBATCH -J {_job_id}
#SBATCH -o {jg.p_job_stdout_filepath}
#SBATCH -e {jg.p_job_stderr_filepath}
#SBATCH -D {jg.p_job_dirpath}
#SBATCH --nodes={p.num_nodes}
#SBATCH --ntasks-per-node={p.num_ranks_per_node}
#SBATCH --mail-type=end 
#SBATCH --mail-user={config_per_user["user_email"]}
#SBATCH --export=NONE 
#SBATCH --time="""+time_str+"""
#SBATCH --cluster=inter
#SBATCH --partition=teramem_inter
#SBATCH --cpus-per-task=94
"""



    content += "\n"
    content += "module load slurm_setup\n"
 

    if False:
    	if p.force_turbo_off:
    		content += """# Try to avoid slowing down some CPUs
#SBATCH --ear=off
"""

    content += """
source /etc/profile.d/modules.sh

module load intel-oneapi-compilers intel-mkl intel-mpi

"""

    if jg.compile.threading != "off":
    	content += """
#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
export OMP_NUM_THREADS="""+str(p.omp_num_threads)+"""
"""

    if p.core_oversubscription:
    	raise Exception("Not supported with this script!")

    if p.core_affinity != None:
    	
    	content += "\necho \"Affnity: "+str(p.core_affinity)+"\"\n"
    	if p.core_affinity == "compact":
    		content += "\nexport OMP_PROC_BIND=close\n"
    	elif p.core_affinity == "scatter":
    		content += "\nexport OMP_PROC_BIND=spread\n"
    	else:
    		raise Exception("Affinity ""+str(p.core_affinity)+"" not supported")

    return content






def jobscript_get_exec_prefix(jg : JobGeneration):
    """
    Prefix before executable

    Returns
    -------
    string
    	multiline text for scripts
    """

    content = ""
    content += jg.runtime.get_jobscript_plan_exec_prefix(jg.compile, jg.runtime)

    return content



def jobscript_get_exec_command(jg : JobGeneration):
    """
    Prefix to executable command

    Returns
    -------
    string
    	multiline text for scripts
    """

    p = jg.parallelization

    mpiexec = ""

    #
    # Only use MPI exec if we are allowed to do so
    # We shouldn"t use mpiexec for validation scripts
    #
    #if not p.mpiexec_disabled:
    #	mpiexec = "mpiexec -n "+str(p.num_ranks)+" --perhost "+str(p.num_ranks_per_node)


    sweet_ld_library_path = os.getenv("MULE_LD_LIBRARY_PATH")
    if sweet_ld_library_path == None:
        raise Exception("Environment variable MULE_LD_LIBRARY_PATH not found!")


    content = f"""

# Output MPI version
echo "**************************************************"
echo "MPI Information"
echo "**************************************************"
echo "mpiexec --version"
mpiexec --version 2>&1
echo "**************************************************"


# List loaded modules
echo "**************************************************"
echo "Loaded modules"
echo "**************************************************"
echo "module list"
module list 2>&1
echo "**************************************************"

# Make sure that MULE library path is really known
export LD_LIBRARY_PATH="{sweet_ld_library_path}:$LD_LIBRARY_PATH"
# mpiexec ... would be here without a line break
EXEC="{jg.compile.getProgramPath()}"
PARAMS="{jg.runtime.getRuntimeOptions()}"
echo "EXEC: ${{EXEC}}"
echo "PARAMS: ${{PARAMS}}"
echo "SLURM_NTASKS: ${{SLURM_NTASKS}}"

{mpiexec} $EXEC $PARAMS || exit 1
"""

    return content



def jobscript_get_exec_suffix(jg : JobGeneration):
    """
    Suffix before executable

    Returns
    -------
    string
    	multiline text for scripts
    """

    content = ""
    content += jg.runtime.get_jobscript_plan_exec_suffix(jg.compile, jg.runtime)
    return content



def jobscript_get_footer(jg : JobGeneration):
    """
    Footer at very end of job script

    Returns
    -------
    string
    	multiline text for scripts
    """

    content = ""
    return content



def jobscript_get_compile_command(jg : JobGeneration):
    """
    Compile command(s)

    This is separated here to put it either
    * into the job script (handy for workstations)
    or
    * into a separate compile file (handy for clusters)

    Returns
    -------
    string
    	multiline text with compile command to generate executable
    """

    content = f"""

SCONS="scons {jg.compile.getSConsParams()} -j 4"
echo "$SCONS"
$SCONS || exit 1
"""

    return content

