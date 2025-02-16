#! /usr/bin/env python3
# Tests space parallelization for ETDSDC on galewsky medium-range (2 days)
import sys
from itertools import product

from mule.JobMule import *
from mule.utils import exec_program
from mule.InfoError import *
from mule.parHelper import *
from mule.etdsdc import getETDSDCSetup

jg = JobGeneration()


"""
Compile parameters
"""
params_compile_sweet_mpi = 'enable'
params_compile_threading = 'omp'
params_compile_thread_parallel_sum = 'enable'

jg.compile.program = 'programs/PDE_SWESphere2D'
# turn off compilation
jg.compilecommand_in_jobscript = False
# !!!!!!!!! time runs
jg.compile.benchmark_timings = 'enable'

jg.compile.cart2d_spectral_space = 'disable'
jg.compile.cart2d_spectral_dealiasing = 'disable'
jg.compile.sphere2d_spectral_space = 'enable'
jg.compile.sphere2d_spectral_dealiasing = 'enable'

jg.unique_id_filter = ['runtime.simparams', 'parallelization', 'runtime.benchmark']


"""
Runtime parameters
"""
jg.runtime.benchmark_name = 'galewsky'

jg.runtime.space_res_spectral = 256
jg.runtime.space_res_grid = None
# NOTE: galewsky does not have analytical solution
jg.runtime.compute_errors = 0

jg.runtime.instability_checks = 0
jg.runtime.verbosity = 2

# Set the number of SH transformation threads
jg.runtime.sh_setup_num_threads = None

"""
Parallelization parameters
"""

# Update TIME parallelization
ptime = JobParallelizationDimOptions('time')
ptime.num_cores_per_rank = 1
ptime.num_threads_per_rank = 1
ptime.num_ranks = 1

pspace = JobParallelizationDimOptions('space')
pspace.num_cores_per_rank = 1
num_threads_per_rank = [1]
pspace.num_ranks = 1

ts_methods = ["ETDSDC(lg,ADDT(lc,n))"]
# test order 3: 
#   optimal (nNodes, nIter)
test_params = [(1, 2)]

#
# Reference solution - we don't really need it here...
#
ref_ts_method = "ERK(ln,order=4)"

jg.runtime.rexi_method = None
jg.runtime.timestepping_method = ref_ts_method
jg.runtime.timestepping_order = 4
jg.runtime.timestepping_order2 = 4
jg.runtime.timestep_size = 15
total_simtime = 1000 * 15 # 1000 steps
jg.runtime.max_simulation_time = total_simtime
jg.runtime.output_timestep_size = jg.runtime.max_simulation_time

jg.reference_job = True
jg.gen_jobscript_directory()
jg.reference_job = False

# Use this one as the reference solution!
jg.reference_job_unique_id = jg.job_unique_id

if __name__ == "__main__":
    for tsm in ts_methods:
        jg.runtime.timestepping_method = tsm
        # generate params for each order
        for (nodes, iters) in test_params:
            SDCparams= getETDSDCSetup(nNodes=nodes, nIter=iters, nodeDistr="LEGENDRE")

            jg.runtime.timestepping_method = tsm
            jg.runtime.paramsSDC = SDCparams
            jg.runtime.init_phase = True
            jg.runtime.nodeType = "LEGENDRE"
            jg.runtime.nNodes = nodes
            jg.runtime.nIter = iters
            jg.runtime.sdcOrder = iters+1
            jg.runtime.idString = tsm +'_'+SDCparams['idString']
            jg.runtime.init_phase = False

            for pspace.num_threads_per_rank in num_threads_per_rank:
                # hack for different job bench folders
                jg.runtime.timestep_size = 60 + pspace.num_threads_per_rank
                # run 1000 steps of each
                jg.runtime.max_simulation_time = 1000 * jg.runtime.timestep_size
                jg.runtime.output_timestep_size = jg.runtime.max_simulation_time
                jg.runtime.sh_setup_num_threads = pspace.num_threads_per_rank
                jg.setup_parallelization([pspace, ptime])
                jg.gen_jobscript_directory()

