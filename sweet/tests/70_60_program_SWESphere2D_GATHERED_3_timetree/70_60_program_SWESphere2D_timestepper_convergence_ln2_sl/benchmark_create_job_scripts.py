#! /usr/bin/env python3

import sys
import copy
from mule.parHelper import *

from mule.JobMule import *
jg = JobGeneration()

from mule.rexi.REXICoefficients import *
from mule.rexi.trexi.TREXI import *
from mule.rexi.cirexi.CIREXI import *
from mule.rexi.brexi.BREXI import *



#
# Run simulation on cart2d or sphere2d
#
jg.compile.program = 'programs/PDE_SWESphere2D'

jg.unique_id_filter = ['runtime.simparams', 'parallelization', 'runtime.benchmark', 'runtime.rexi_params']


# Verbosity mode
jg.runtime.verbosity = 2

#
# Mode and Grid resolution
#
jg.runtime.space_res_spectral = 64
jg.runtime.space_res_grid = None

#jg.runtime.benchmark_name = "gaussian_bumps_test_cases"
jg.runtime.benchmark_name = "galewsky"


#
# Compute error
#
jg.runtime.compute_errors = 0

#
# Preallocate the REXI matrices
#
jg.runtime.rexi_sphere2d_preallocation = 1

#
# Threading across all REXI terms
#
jg.compile.rexi_thread_parallel_sum = "enable"



jg.runtime.f_sphere2d = 0

jg.runtime.viscosity = 0.0



#ref_ts_size = 8
ref_ts_method = f"ERK(ln,order=4)"

tso = 2
ts_methods = [
            f"SETTLS(lg,na(sl_order=2),ADDT(lc,nr),order=2)",
            ###f"SLETDRK(EXP(lg),na(sl_order=2),ADDT(lc,nr),order=1)",
    ]


timestep_size_min = 16
timestep_sizes = [timestep_size_min*(2.0**i) for i in range(0, 5)]

jg.runtime.max_simulation_time = timestep_size_min*512

#####################################################
#####################################################
#####################################################

jg.runtime.output_timestep_size = jg.runtime.max_simulation_time

jg_base = copy.deepcopy(jg)

jg = None






for tsm in ts_methods:
    for timestep_size in timestep_sizes:

        space_res_spectral = 32 * (timestep_sizes[-1] // timestep_size)

        #
        # Reference solution for particular resolution
        #

        jgref = copy.deepcopy(jg_base)
        jgref.runtime.rexi_method = None
        jgref.runtime.timestepping_method = ref_ts_method
        #jgref.runtime.timestep_size = ref_ts_size
        jgref.runtime.timestep_size = timestep_size//4

        jgref.runtime.space_res_spectral = space_res_spectral

        jgref.reference_job = True
        jgref.gen_jobscript_directory()
        jgref.reference_job = False
        # Use this one as the reference solution!

        #
        # Regular job
        #
        jg = copy.deepcopy(jg_base)

        jg.reference_job_unique_id = jgref.job_unique_id

        jg.runtime.timestepping_method = tsm
        jg.runtime.timestepping_order = tso
        jg.runtime.timestepping_order2 = tso
        jg.runtime.timestep_size = timestep_size

        jg.runtime.space_res_spectral = space_res_spectral

        setupParallelization(jg)
        jg.gen_jobscript_directory()
