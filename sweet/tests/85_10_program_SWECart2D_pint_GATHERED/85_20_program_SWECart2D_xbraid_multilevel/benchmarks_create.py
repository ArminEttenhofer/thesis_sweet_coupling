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

tsm_fine = "SS(IRK(l,order=2),ERK(n,order=2),order=2)"
###tsm_fine = "SLETDRK(EXP(l),na(sl_order=2),nr,order=2)"

#Create main compile/run options
jg = JobGeneration()

#Get Earth parameters (if necessary)
earth = EarthMKSDimensions()


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
##jg.runtime.timestepping_order = orders[tsm_fine]
##jg.runtime.timestepping_order2 = orders[tsm_fine]
jg.runtime.space_res_spectral = 32


## Reference job
jg.reference_job = True

jg.compile.program = "programs/PDE_SWECart2D"
jg.compile.mode = "debug"
jg.compile.sweet_mpi = "disable"
jg.compile.xbraid_scalar = 'disable'

jg.compile.xbraid = "none";
jg.runtime.xbraid_enabled = 0;

jg.compile.parareal = "none";
jg.runtime.parareal_enabled = 0;

jg.runtime.timestepping_method = tsm_fine
###jg.runtime.timestepping_order = orders[tsm_fine]
###jg.runtime.timestepping_order2 = orders[tsm_fine]

jg.gen_jobscript_directory();


## MGRIT jobs

jg.reference_job = False
jg.reference_job_unique_id = jg.job_unique_id
ref_path = jg.p_job_dirpath
print("REF PATH ", ref_path)

jg.compile.program = "programs/xbraid_PDE_SWECart2D"
jg.compile.mode = "debug"
jg.compile.sweet_mpi = "enable"
jg.compile.xbraid_cart2d = 'enable'
jg.compile.xbraid_cart2d_swe = 'enable'

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
jg.runtime.xbraid_max_iter = 10
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
jg.runtime.xbraid_timestepping_order = 2
jg.runtime.xbraid_timestepping_order2 = 2
jg.runtime.xbraid_verbosity = 0;
jg.runtime.xbraid_load_ref_csv_files = 0;
jg.runtime.xbraid_path_ref_csv_files = "";
jg.runtime.xbraid_load_fine_csv_files = 0;
jg.runtime.xbraid_path_fine_csv_files = "";
jg.runtime.xbraid_store_iterations = 0;
jg.runtime.xbraid_spatial_coarsening = 0;
jg.runtime.xbraid_pt = 1;

jg.runtime.xbraid_max_levels = 1
jg.runtime.xbraid_store_iterations = 1;

tsms_fine = [tsm_fine]
##tsms_coarse = ["l_rexi_n_erk", "l_rexi_n_etdrk", "l_cn_na_sl_nd_settls"]
###tsms_coarse = [tsm_fine,"SLETDRK(EXP(l),na(sl_order=2),nr,order=2)"]
###tsms_coarse = ["SLETDRK(EXP(l),na(sl_order=2),nr,order=2)"]
###tsms_coarse = ["SS(IRK(l,order=2),ERK(n,order=2),order=2)", "SLETDRK(EXP(l),na(sl_order=2),nr,order=2,version=default)"]
##tsms_coarse = ["SS(IRK(l,order=2),ERK(n,order=2),order=2)", "SETTLS(l,na(sl_order=2),nr,order=2)"]
tsms_coarse = ["SS(IRK(l,order=2),ERK(n,order=2),order=2)"]
jg.runtime.xbraid_store_iterations = 1;
jg.runtime.xbraid_access_level = 2;

cfactors = [2, 4];
nbs_levels = [2, 4];
nb_pts = [1, 2];
spatial_coarsenings = [1]
##nb_pts = [1, 2, 4];

for tsm_fine in tsms_fine:
    for tsm_coarse in tsms_coarse:
        for online_error in range(2):
            for cfactor in cfactors:
                for nb_levels in nbs_levels:
                    for pt in nb_pts:
                        for spatial_coarsening in spatial_coarsenings:

                            jg.runtime.xbraid_timestepping_method = tsm_fine + ";" + tsm_coarse;
                            ####jg.runtime.xbraid_timestepping_order = "2";
                            ####jg.runtime.xbraid_timestepping_order2 = "2";

                            if online_error:
                                jg.runtime.xbraid_store_iterations = 0;
                                jg.runtime.xbraid_load_fine_csv_files = 1;
                                jg.runtime.xbraid_path_fine_csv_files = ref_path;
                            else:
                                jg.runtime.xbraid_store_iterations = 1;
                                jg.runtime.xbraid_load_fine_csv_files = 0;

                            jg.runtime.xbraid_cfactor = cfactor;
                            jg.runtime.xbraid_max_levels = nb_levels;
                            jg.runtime.xbraid_spatial_coarsening = spatial_coarsening;
                            jg.runtime.xbraid_pt = pt;

                            params_pspace_num_cores_per_rank = [jg.platform_resources.num_cores_per_socket]
                            params_pspace_num_threads_per_rank = [jg.platform_resources.num_cores_per_socket]
                            params_ptime_num_cores_per_rank = [1]

                            # Update TIME parallelization
                            ptime = JobParallelizationDimOptions('time')
                            ptime.num_cores_per_rank = 1
                            ptime.num_threads_per_rank = 1 #pspace.num_cores_per_rank
                            ptime.num_ranks = pt

                            pspace = JobParallelizationDimOptions('space')
                            pspace.num_cores_per_rank = 1
                            pspace.num_threads_per_rank = params_pspace_num_cores_per_rank[-1]
                            pspace.num_ranks = 1

                            # Setup parallelization
                            jg.setup_parallelization([pspace, ptime], override_insufficient_resources=True)


                            jg.gen_jobscript_directory();
