#! /usr/bin/env python3

import os
import sys
import math
import copy

from itertools import product

from mule.JobGeneration import *
from mule.JobMule import *



jg = JobGeneration()

#
# Compile options
#
jg.compile.program = "programs/PDE_SWESphere2D"

jg.compile.threading = "omp"

jg.unique_id_filter += ["runtime.benchmark", "runtime"]


#
# Runtime options
#
#jg.runtime.space_res_spectral = 256
jg.runtime.space_res_spectral = 128

#jg.runtime.max_simulation_time = 60*60*24*8    # 8 days
jg.runtime.max_simulation_time = 60*60*1    # 1 hour
jg.runtime.output_timestep_size = -1
jg.runtime.benchmark_name = "galewsky"
jg.runtime.reuse_plans = "load" # Try to load plans, if not, create them

# Saves us some time in case of unstable simulations
jg.runtime.instability_checks = 0

jg.runtime.timestepping_method = "ERK(ln,order=4)"

jg.runtime.timestep_size = 30

jg.parallelization.max_wallclock_seconds = 1*60*60

jg.compilecommand_in_jobscript = False

jgbase = copy.deepcopy(jg)
jg = None


compile_commands_accum = []


def setupParallelization(jg, num_threads = None):

    if num_threads is None:
        num_threads = jgbase.platform_resources.num_cores_per_node

    pspace = JobParallelizationDimOptions()
    pspace.num_cores_per_rank = jgbase.platform_resources.num_cores_per_node
    pspace.num_threads_per_rank = num_threads
    pspace.num_ranks = 1

    jg.runtime.sh_setup_num_threads = num_threads

    jg.setup_parallelization([pspace])
    jg.parallelization.core_oversubscription = False

    jg.parallelization.core_affinity = "compact"



#
# Generate some solution which produces output
#
if True:
    jgoutput = copy.deepcopy(jgbase)

    setupParallelization(jgoutput)

    jgoutput.runtime.output_timestep_size = jgoutput.runtime.max_simulation_time

    jgoutput.gen_jobscript_directory("job_output_"+jgoutput.getUniqueID())

    for i in jgoutput.get_compilecommands_accum():
        if i not in compile_commands_accum:
            compile_commands_accum.append(i)




#
# Create job scripts
#
if True:

    for i in range(1, jgbase.platform_resources.num_cores_per_node):
        jg = copy.deepcopy(jgbase)

        setupParallelization(jg, i)
        
        jg.gen_jobscript_directory("job_bench_"+jg.getUniqueID())

        for i in jg.get_compilecommands_accum():
            if i not in compile_commands_accum:
                compile_commands_accum.append(i)

jg.write_compilecommands(content_accum=compile_commands_accum)

