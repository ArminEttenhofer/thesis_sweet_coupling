#! /usr/bin/env python3

from mule.JobGeneration import *
jg = JobGeneration()

jg.runtime.floating_point_output_digits = 16

jg.compile.program = 'programs/PDE_SWESphere2D'
jg.compile.mode = 'release'
#jg.compile.mode = 'debug'

jg.compile.cart2d_spectral_space = 'disable'
jg.compile.cart2d_spectral_dealiasing = 'disable'
jg.compile.sphere2d_spectral_space = 'enable'
jg.compile.sphere2d_spectral_dealiasing = 'enable'

# Request csv output for comparisons
jg.runtime.output_file_mode = 'csv'

jg.runtime.space_res_spectral = "66,66"
jg.runtime.space_res_grid = "196,100"
jg.runtime.instability_checks = 0
jg.runtime.benchmark_name = 'galewsky'

for analytical in [0, 1]:

    #jg.compile.threading = 'omp'
    jg.compile.threading = 'off'

    jg.unique_id_filter = ['runtime', 'compile', 'parallelization']

    #jg.runtime.gravitation= 1
    #jg.runtime.sphere2d_rotating_coriolis_omega = 1
    #jg.runtime.h0 = 1
    #jg.runtime.cart2d_domain_size = 1

    jg.runtime.verbosity = 10

    if False:
            jg.runtime.timestep_size = 60
            jg.runtime.max_simulation_time = 60*10
            jg.runtime.output_timestep_size = 60
    else:
            jg.runtime.timestep_size = 10
            jg.runtime.max_simulation_time = 60*60*24
            jg.runtime.output_timestep_size = 60*60

    jg.runtime.output_time_scale_inv = 60*60

    jg.runtime.output_filename = "output_%s_t%020.8f.csv"

    if analytical:
        jg.runtime.benchmark_galewsky_geostropic_setup = "analytical"
        file_ext = "_analytical_1"
    else:
        jg.runtime.benchmark_galewsky_geostropic_setup = "numerical"
        file_ext = "_analytical_0"

    jg.runtime.timestepping_method = 'ERK(ln,o=1)'

    jg.gen_jobscript_directory('job_bench_sweet'+file_ext)

