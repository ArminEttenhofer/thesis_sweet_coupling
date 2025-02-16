#! /usr/bin/env python3

import os
import sys
import math
import copy

from itertools import product

from mule.JobGeneration import *

jg = JobGeneration()

#
# Compile options
#
jg.compile.program = 'programs/PDE_SWESphere2D'

#
# Runtime options
#
#jg.runtime.space_res_spectral = 256
jg.runtime.space_res_spectral = 64

#jg.runtime.max_simulation_time = 60*60*24*8    # 8 days
jg.runtime.max_simulation_time = 60*60*24*1    # 1 day
jg.runtime.output_timestep_size = jg.runtime.max_simulation_time
jg.runtime.benchmark_name = "galewsky"

# Saves us some time in case of unstable simulations
jg.runtime.instability_checks = 1

jgbase = copy.deepcopy(jg)
jg = None


#
# Which benchmarks to run
#

# Timestep sizes
#params_timestep_sizes_ = [15/8, 15/4, 15/2, 15, 30, 60, 120, 180, 360, 720]
params_timestep_sizes_ = [30, 60, 120, 180, 360, 720]

# Time stepping methods
ts_methods = [
    'ERK(ln,order=2)',
    'SS(IRK(lg,order=2),ERK(ADDT(lc,n),order=2),order=2)',
    'SS(IRK(l,order=2),ERK(n,order=2),order=2)'
]



#
# Reference benchark
#

# Ts parameter
ref_ts_method = 'ERK(ln,order=4)'   # Used as reference solution

# Pick the smallest time step size for the reference time step size
timestep_size_reference = params_timestep_sizes_[0]




#
# SHTNS plan generation scripts
#
if True:
    jgplan = copy.deepcopy(jgbase)

    # We require the plans to be recreated
    # This allows better comparison of different runs since the plans could otherwise change between runs
    jgplan.runtime.reuse_plans = "save"

    # Just some dummy time step size to get things going

    jgplan.runtime.timestepping_method = ref_ts_method

    jgplan.runtime.timestep_size = params_timestep_sizes_[0]

    # Set simtime to 0
    jgplan.runtime.max_simulation_time = 0

    # Skip postprocessing for plan generation
    jgplan.skip_postprocessing = True

    # No output
    jgplan.runtime.output_timestep_size = -1
    jgplan.runtime.output_filename = "-"

    # We need to use a special prefix for this job
    jobdir = 'job_plan'+jgplan.getUniqueID()

    # Generate jobs
    jgplan.gen_jobscript_directory(jobdir)

    # Now we can reuse these plans
    jgplan.runtime.reuse_plans = "require_load"





#
# Generate reference solution
#
if True:
    jgref = copy.deepcopy(jgbase)

    jgref.runtime.timestep_size  = timestep_size_reference

    jgref.runtime.timestepping_method = ref_ts_method

    # Set this to true to say that this is one of the reference jobs
    jgref.reference_job = True

    jgref.gen_jobscript_directory('job_benchref_'+jgref.getUniqueID())






#
# Create job scripts
#
for tsm in ts_methods:
    for ts in params_timestep_sizes_:

        jg = copy.deepcopy(jgbase)

        jg.runtime.timestep_size = ts
        jg.runtime.timestepping_method = tsm

        # We now reuse the unique job ID of the reference solution to tell the other jobs about their reference solution!
        jg.reference_job_unique_id = jgref.job_unique_id

        jg.gen_jobscript_directory('job_bench_'+jg.getUniqueID())

