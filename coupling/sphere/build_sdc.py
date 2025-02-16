#! /usr/bin/env python3
import copy

from mule.JobGeneration import JobGeneration
from mule.JobParallelizationDimOptions import JobParallelizationDimOptions
from mule.sdc import getSDCSetup

jgbase = JobGeneration()
verbose = True

jgbase.compile.mode = "release"
jgbase.compile.gui = "disable"

#
# Mode and Physical resolution
#
jgbase.runtime.space_res_spectral = 2501
jgbase.runtime.space_res_grid = -1
nProcSpace = 8

jgbase.parallelization.core_oversubscription = False
jgbase.parallelization.core_affinity = "compact"

jgbase.compile.threading = "omp"
jgbase.compile.rexi_thread_parallel_sum = "disable"

gen_reference_solution = False
jgbase.runtime.benchmark_name = "column"

jgbase.runtime.max_simulation_time = 1000

jgbase.runtime.output_timestep_size = 5
jgbase.runtime.output_file_mode = "bin"

base_timestep_size = 0.5
sdcParameters = dict(
    nNodes=4,
    nIter=3,
    nodeType='RADAU-RIGHT', 
    qDeltaImplicit='MIN-SR-S', 
    qDeltaExplicit='PIC', 
    qDeltaInitial='BEPAR',
    preSweep="QDELTA",
    postSweep="LASTNODE",
    nodeDistr='LEGENDRE',
)
jgbase.runtime.timestepping_method = "SDC-FP(i=lg,ADDT(lc,n))"
# jgbase.runtime.timestepping_method = "ERK(ln,o=4)"
# jgbase.runtime.viscosity = 0.5

params_pspace_num_cores_per_rank = [jgbase.platform_resources.num_cores_per_socket]
params_pspace_num_threads_per_rank = [jgbase.platform_resources.num_cores_per_socket]

jgbase.parallelization.max_wallclock_seconds = 60*60*1  # Limit to one hour


# Resolution specific multiplier
resmul = 128/jgbase.runtime.space_res_spectral

jgbase.unique_id_filter = [
    "compile",
    "runtime.max_simulation_time",
    "runtime.benchmark",
    "compile_cart2d",
]
jgbase.compilecommand_in_jobscript = True

compile_commands_accum = []
def compileCommandsAccum(jg):
    for i in jg.get_compilecommands_accum():
        if i not in compile_commands_accum:
            compile_commands_accum.append(i)

#
# Run simulation on cart2d or sphere
#
jgbase.compile.program = "programs/PDE_SWESphere2D"

jgbase.compile.cart2d_spectral_space = "enable"
jgbase.compile.cart2d_spectral_dealiasing = "disable"
jgbase.compile.sphere2d_spectral_space = "enable"
jgbase.compile.sphere2d_spectral_dealiasing = "enable"

jgbase.runtime.verbosity = 0
jgbase.runtime.instability_checks = 0

def setupParallelizationSDC(jg, ranks_in_time=None, threads_in_time=None):


    pspace = JobParallelizationDimOptions("space")
    pspace.num_cores_per_rank = nProcSpace # jgbase.platform_resources.num_cores_per_socket
    pspace.num_threads_per_rank = nProcSpace # jgbase.platform_resources.num_cores_per_socket
    pspace.num_ranks = 1
    #pspace.threading = "omp"

    jg.runtime.sh_setup_num_threads = pspace.num_threads_per_rank
    jg.parallelization.withNestedOpenMPNumThreads = False

    if ranks_in_time is not None or threads_in_time is not None:
        # Update TIME parallelization
        ptime = JobParallelizationDimOptions("time")
        ptime.num_cores_per_rank = threads_in_time
        ptime.num_threads_per_rank = threads_in_time
        ptime.num_ranks = ranks_in_time
        #ptime.threading = "omp"

        ptime.print()

        jg.setup_parallelization([ptime, pspace])

    else:
        jg.setup_parallelization([pspace])

    if verbose:
        print("Parallelization output:")
        jg.parallelization.print()


jgsdc = copy.deepcopy(jgbase)


# We first setup all parameters
jgsdc.runtime.timestep_size = base_timestep_size

paramsSDC = getSDCSetup(**sdcParameters)
jgsdc.runtime.sdc_params = paramsSDC

jgsdc.runtime.sdc_parallel = 0
setupParallelizationSDC(jgsdc, 1, 1)

# Ask for jobID (including all parameters)
jobdir = "job_bench_"+jgsdc.getUniqueID()

# Set the SDC file
jgsdc.runtime.sdc_file = "params_SDC.sweet"

jgsdc.gen_jobscript_directory(jobdir)

# Write out SDC coefficients after job directory has been created
paramsSDC.writeToFile(jgsdc.job_dirpath+"/"+jgsdc.runtime.sdc_file)

compileCommandsAccum(jgsdc)




jgsdc.write_compilecommands(content_accum=compile_commands_accum)


