#! /usr/bin/env python3
#
#  Create series of job to be run with sweet
#
#  Pedro Peixoto <pedrosp@ime.usp.br>
#  modified from Martin Schreiber initial job_create
#
#-------------------------------------------------------

import os
import sys
import stat
import math
from glob import glob

#Classes containing sweet compile/run basic option
from mule.JobGeneration import *
from mule.SWEETRuntimeParametersScenarios import *
from mule.JobParallelization import *
from mule.JobParallelizationDimOptions import *

orders = {};
orders["l_irk_n_erk"] = 1;
orders["l_irk"] = 1;
orders["l_erk_n_erk"] = 2;
orders["l_erk"] = 2;
orders["l_cn_n_erk"] = 2;
orders["l_cn"] = 2;
orders["l_rexi_n_etdrk"] = 2;
orders["l_rexi_n_erk"] = 2;
orders["l_rexi"] = 2;
orders["ln_erk"] = 2;
orders["l_direct"] = 2;
orders["l_rexi_na_sl_nd_settls"] = 2;
orders["l_rexi_na_sl_nd_etdrk"] = 2;
orders["l_cn_na_sl_nd_settls"] = 2;

tsm_fine = "ln_erk"
tsm_coarse = "ln_erk"

tsm_fine = "SS(IRK(l,order=2),ERK(n,order=2),order=2)"
tsm_coarse = "SS(IRK(l,order=2),ERK(n,order=2),order=2)"

#Create main compile/run options
jg = JobGeneration()

#Get Earth parameters (if necessary)
earth = EarthMKSDimensions()

#
# Run simulation on cart2d or sphere2d
#
#Basic cart2d options
jg.compile.program = "programs/xbraid_PDE_SWECart2D"
jg.compile.mode = "debug"
###<<<<<<< HEAD:tests/87_10_program_swe_sphere_pint_GATHERED/87_20_program_swe_sphere_xbraid_wrapped_functions/benchmarks_create.py
########jg.compile.sweet_mpi = "enable"
########jg.compile.xbraid_sphere = 'enable'
###jg.compile.fortran_source = "enable"
###jg.compile.lapack = "enable"
###=======
jg.compile.sweet_mpi = "enable"
jg.compile.xbraid_cart2d = 'enable'
jg.compile.xbraid_cart2d_swe = 'enable'
###>>>>>>> main:tests/85_10_program_SWECart2D_pint_GATHERED/85_20_program_SWECart2D_xbraid_wrapped_functions/benchmarks_create.py

jg.runtime.output_file_mode = "csv"

# Verbosity mode
jg.runtime.verbosity = 3

jg.compile.cart2d_spectral_space = "enable";
jg.compile.cart2d_spectral_dealiasing = "enable";

#
# Benchmark ID
#
#jg.runtime.bench_id = 1
jg.runtime.benchmark_name = "unstablejet"

#
# Compute error or difference to initial data
#
jg.runtime.compute_errors = 1

jg = DisableGUI(jg)

#
# REXI method
###jg.runtime.rexi_method = 'direct'
#jg.runtime.rexi_use_direct_solution = 1

jg = RuntimeSWECart2DEarthParam(jg)
#jg = RuntimeSWENonDimParam(jg)

jg.runtime.viscosity = 0.0

# Deactivate threading
jg.compile.threading = "omp"

#
# Time, Mode and Grid resolution
#
jg.runtime.max_simulation_time = 100.
jg.runtime.output_timestep_size = 25.
timestep_size_reference = 2.5
timestep_size_fine = 5.
jg.runtime.timestep_size = timestep_size_fine
jg.runtime.timestepping_method = tsm_fine
if "(" not in tsm_fine:
    jg.runtime.timestepping_order = orders[tsm_fine]
    jg.runtime.timestepping_order2 = orders[tsm_fine]
jg.runtime.space_res_spectral = 32

jg.compile.xbraid = "mpi";
jg.runtime.xbraid_enabled = 1;
jg.runtime.xbraid_max_levels = 3
jg.runtime.xbraid_skip = 1
jg.runtime.xbraid_min_coarse = 2
jg.runtime.xbraid_nrelax = 1
jg.runtime.xbraid_nrelax0 = -1
jg.runtime.xbraid_tol = 1e-14
jg.runtime.xbraid_tnorm = 2
jg.runtime.xbraid_cfactor = 2
jg.runtime.xbraid_cfactor0 = -1
jg.runtime.xbraid_max_iter = 100
jg.runtime.xbraid_fmg = 0
jg.runtime.xbraid_res = 0
jg.runtime.xbraid_storage = 0
jg.runtime.xbraid_print_level = 2
jg.runtime.xbraid_access_level = 1
jg.runtime.xbraid_run_wrapper_tests = 0
jg.runtime.xbraid_fullrnorm = 2
jg.runtime.xbraid_use_seq_soln = 0
jg.runtime.xbraid_use_rand = 1
jg.runtime.xbraid_timestepping_method = tsm_fine
if "(" not in tsm_fine:
    jg.runtime.xbraid_timestepping_order = orders[tsm_fine]
    jg.runtime.xbraid_timestepping_order2 = orders[tsm_fine]
jg.runtime.xbraid_verbosity = 0;
jg.runtime.xbraid_load_ref_csv_files = 0;
jg.runtime.xbraid_path_ref_csv_files = "";
jg.runtime.xbraid_load_fine_csv_files = 0;
jg.runtime.xbraid_path_fine_csv_files = "";
jg.runtime.xbraid_store_iterations = 0;
jg.runtime.xbraid_spatial_coarsening = 0;
jg.runtime.xbraid_pt = 1;

jg.runtime.xbraid_run_wrapper_tests = 1

jg.gen_jobscript_directory();
