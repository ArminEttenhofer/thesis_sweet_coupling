#! /usr/bin/env python3
# Tests convergence rates for ETDSDC on galewsky medium-range (1 day)
import sys
from itertools import product

from mule.JobMule import *
from mule.utils import exec_program
from mule.InfoError import *
from mule.parHelper import *
from mule.etdsdc import getETDSDCSetup
from mule.sdc import getSDCSetup

jg = JobGeneration()


"""
Compile parameters
"""
params_compile_sweet_mpi = 'enable'
params_compile_threading = 'omp'
params_compile_thread_parallel_sum = 'enable'

jg.compile.program = 'programs/PDE_SWESphere2D'
# turn off compilation: coolmuc2
jg.compilecommand_in_jobscript = False
# !!!!!!! time runs
jg.compile.benchmark_timings = 'enable'

jg.compile.cart2d_spectral_space = 'disable'
jg.compile.cart2d_spectral_dealiasing = 'disable'
jg.compile.sphere2d_spectral_space = 'enable'
jg.compile.sphere2d_spectral_dealiasing = 'enable'

jg.unique_id_filter = ['runtime.simparams', 'parallelization', 'runtime.benchmark']

"""
Runtime parameters
"""
params_runtime_timestep_sizes = [30., 60., 120., 240., 360., 480., 600., 960., 1920.] # [s]
total_simtime = 60 * 60 * 24 * 15 # 15 days

jg.runtime.benchmark_name = 'williamson5' # flow over an isolated mountain

jg.runtime.space_res_spectral = 128
jg.runtime.space_res_grid = None
# NOTE: galewsky does not have analytical solution
jg.runtime.compute_errors = 0
jg.runtime.max_simulation_time = total_simtime
# shortcut instable simulations
jg.runtime.instability_checks = 0
jg.runtime.verbosity = 2

# output results after end
jg.runtime.output_timestep_size = jg.runtime.max_simulation_time

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
pspace.num_threads_per_rank = 56 # see strong scaling
pspace.num_ranks = 1

jg.setup_parallelization([pspace, ptime])

ts_methods = ["SDC_Classic(lg,ADDT(lb,lc,n))", "ETDSDC(lg,ADDT(lb,lc,n))"]
# test orders
test_orders = [2, 3, 4]

#
# Reference solution
#
ref_ts_method = "ERK(ln,order=4)"

jg.runtime.rexi_method = None
jg.runtime.timestepping_method = ref_ts_method
jg.runtime.timestepping_order = 4
jg.runtime.timestepping_order2 = 4
jg.runtime.timestep_size = 15

jg.reference_job = True
jg.gen_jobscript_directory()
jg.reference_job = False

# Use this one as the reference solution!
jg.reference_job_unique_id = jg.job_unique_id

if __name__ == "__main__":
    for tsm in ts_methods:
        
        # generate params for each order
        if (tsm.find("ETDSDC") != -1):
            jg.runtime.timestepping_method = tsm
            for order in test_orders:
                ETDSDCparams= getETDSDCSetup(nNodes=max(order-2, 0), nIter=max(order-1, 0), nodeDistr="LEGENDRE")

                jg.runtime.paramsSDC = ETDSDCparams
                jg.runtime.init_phase = True
                jg.runtime.idString = tsm +'_'+ETDSDCparams['idString']
                jg.runtime.nodeType = "LEGENDRE"
                jg.runtime.nNodes = max(order-2, 0)
                jg.runtime.nIter = max(order-1, 0)
                jg.runtime.sdcOrder = order
                jg.runtime.init_phase = False
                
                for timestep in params_runtime_timestep_sizes:
                    jg.runtime.timestep_size = timestep
                    jg.gen_jobscript_directory()
        
        elif (tsm.find("Classic") != -1):
            jg.runtime.timestepping_method = tsm
            for order in test_orders:
                SDCparams= getSDCSetup(nNodes=max(order, 2), nIter=max(order-1, 0), nodeDistr="LEGENDRE", nodeType="LOBATTO")

                jg.runtime.paramsSDC = SDCparams
                jg.runtime.init_phase = True
                jg.runtime.idString = tsm +'_'+SDCparams['idString']
                jg.runtime.nodeType = "LEGENDRE"
                jg.runtime.nNodes = max(order, 2)
                jg.runtime.nIter = max(order-1, 0)
                jg.runtime.sdcOrder = order
                jg.runtime.init_phase = False

                for timestep in params_runtime_timestep_sizes:
                    jg.runtime.timestep_size = timestep
                    jg.gen_jobscript_directory()

