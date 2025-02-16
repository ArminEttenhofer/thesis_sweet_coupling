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
jg.compile.mode = 'debug'

jg.unique_id_filter = ['runtime.simparams', 'parallelization', 'runtime.benchmark', 'runtime.rexi_params']


# Verbosity mode
jg.runtime.verbosity = 2

#
# Mode and Grid resolution
#
jg.runtime.space_res_spectral = 64
jg.runtime.space_res_grid = None

jg.runtime.benchmark_name = "gaussian_bumps_test_cases"
#jg.runtime.benchmark_name = "galewsky"

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



ref_ts_size = 2
ref_ts_order = 4
o=f"order={ref_ts_order}"
ref_ts_method = f"ERK(ln,{o})"
#ref_ts_method = "ln_erk"

ts_order = 2
o = f"order={ts_order}"
ts_methods = [
            #'ln_erk',
            f"ERK(ln,{o})",

            #"l_exp_n_erk_ver0",
            f"SS(REXI(l),ERK(n,{o}),{o})",
            #"l_exp_n_erk_ver1",
            f"SS(ERK(n,{o}),REXI(l),{o})",

            #"lg_exp_lc_n_erk_ver0",
            f"SS(EXP(lg),ERK(ADDT(lc,n),{o}),{o})",
            f"SS(REXI(lg),ERK(ADDT(lc,n),{o}),{o})",
            #"lg_exp_lc_n_erk_ver1",
            f"SS(ERK(ADDT(lc,n),{o}),EXP(lg),{o})",
            f"SS(ERK(ADDT(lc,n),{o}),REXI(lg),{o})",

            #"l_exp_n_etdrk",
            f"ETDRK(REXI(l),n,{o})",

            f"ETDRK(REXI(lg),ADDT(lc,n),{o})",
            f"ETDRK(EXP(lg),ADDT(lc,n),{o})",

            #"lg_exp_na_sl_lc_nr_etdrk_uv",
            f"SLETDRK(EXP(lg),na(sl_order=2),ADDT(lc,nr),{o})",
    ]


timestep_size_min = 32
timestep_sizes = [timestep_size_min*(2.0**i) for i in range(0, 6)]

jg.runtime.max_simulation_time = timestep_size_min*512

#####################################################
#####################################################
#####################################################

jg.runtime.output_timestep_size = jg.runtime.max_simulation_time



#
# Reference solution
#
jg.runtime.rexi_method = None
jg.runtime.timestepping_method = ref_ts_method
jg.runtime.timestepping_order = ref_ts_order
jg.runtime.timestepping_order2 = ref_ts_order
jg.runtime.timestep_size = ref_ts_size

jg.reference_job = True
jg.gen_jobscript_directory()
jg.reference_job = False

# Use this one as the reference solution!
jg.reference_job_unique_id = jg.job_unique_id



#
# Create job scripts
#

jg.runtime.timestepping_order = ts_order
jg.runtime.timestepping_order2 = ts_order


jg_base = copy.deepcopy(jg)

for tsm in ts_methods:
    for timestep_size in timestep_sizes:

        # Create new job description
        jg = copy.deepcopy(jg_base)

        jg.runtime.timestepping_method = tsm
        jg.runtime.timestep_size = timestep_size

        if jg.runtime.max_simulation_time % jg.runtime.timestep_size != 0:
            print("simtime: "+str(jg.runtime.max_simulation_time))
            print("timestep_size: "+str(jg.runtime.timestep_size))
            raise Exception("Invalid time step size (not remainder-less dividable)")

        if 'REXI(' in jg.runtime.timestepping_method:

            if 'rexi' in jg.runtime.timestepping_method.lower():
                # CI REXI via file input
                jg.runtime.rexi_method = 'file'

                #
                # !!! WARNING !!!
                # A timestep size of 32 leads to issues with the REXI method on Gitlab CI.
                # Probably, the solver accuracy is degraded due to the accummulation of the REXI terms
                #
                jg.runtime.timestep_size *= 2

                cirexi = CIREXI()
                coeffs_phi0 = cirexi.setup("phi0", N=16, R=2).toFloat()
                coeffs_phi1 = cirexi.setup("phi1", N=16, R=2).toFloat()
                coeffs_phi2 = cirexi.setup("phi2", N=16, R=2).toFloat()
                coeffs_phi3 = cirexi.setup("phi3", N=16, R=2).toFloat()

                coeffs_phi0.normalize_steady_state()
                coeffs_phi1.normalize_steady_state()
                coeffs_phi2.normalize_steady_state()
                coeffs_phi3.normalize_steady_state()

                jg.runtime.rexi_files_coefficients = [coeffs_phi0, coeffs_phi1, coeffs_phi2, coeffs_phi3]

            else:
                # B REXI via file
                jg.runtime.rexi_method = 'file'

                brexi = BREXI()
                coeffs = brexi.setup(N=8, quadrature_method='gauss_legendre').toFloat()
                jg.runtime.rexi_files_coefficients = [coeffs]


        setupParallelization(jg)
        jg.gen_jobscript_directory()
