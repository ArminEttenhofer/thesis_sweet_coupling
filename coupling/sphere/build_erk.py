#! /usr/bin/env python3

import os
import sys
import stat
import math

from mule.JobMule import *
jg = JobGeneration()

###################################################
# Compilation Settings for reference & tests jobs #
###################################################

jg.compile.program = 'programs/PDE_SWESphere2D'

# run simulation on sphere2d, not cart2d
jg.compile.cart2d_spectral_space = 'disable'
jg.compile.cart2d_spectral_dealiasing = 'disable'
jg.compile.sphere2d_spectral_space = 'enable'
jg.compile.sphere2d_spectral_dealiasing = 'enable'

# enable MPI
jg.compile.sweet_mpi = 'disable'

# other compilation settings
jg.compile.libsph = 'enable'
jg.compile.threading = 'off'
jg.compile.libfft = 'enable'
#jg.compile.quadmath = 'enable'

###############################################
# Runtime Settings for reference & tests jobs #
###############################################

jg.runtime.output_file_mode = 'bin'

# Verbosity mode
jg.runtime.verbosity = 2

# Mode and Grid resolution
jg.runtime.space_res_spectral = 200
jg.runtime.space_res_grid = None

# Benchmark
jg.runtime.benchmark_name = "column"

# Compute error
jg.runtime.compute_errors = 0

jg.runtime.f_sphere2d = 0

jg.unique_id_filter = ['compile', 'parallelization']

jg.runtime.timestepping_method = 'l_erk'
jg.runtime.timestepping_order = 4
jg.runtime.timestepping_order2 = 4

jg.runtime.timestep_size = 1
jg.runtime.max_simulation_time = 1000
jg.runtime.output_timestep_size = 5

#################
# Reference Job #
#################

jg.reference_job = False
jg.gen_jobscript_directory()
