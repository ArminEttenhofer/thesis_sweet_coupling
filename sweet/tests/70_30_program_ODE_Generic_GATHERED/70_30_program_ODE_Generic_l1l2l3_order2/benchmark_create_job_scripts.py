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


jg.unique_id_filter = [
        "compile",
        "parallelization",
        "runtime.benchmark",
    ]

jg.compile.program = "programs/ODE_Generic"
jg.compile.rexi_thread_parallel_sum = "enable"
jg.compile.mode = "release"


jg.runtime.verbosity = 2
jg.runtime.compute_errors = 0
jg.runtime.benchmark_name = "initone"
jg.runtime.ode = "dahlquist"

jg.runtime.max_simulation_time = 100

jg.runtime.output_timestep_size = jg.runtime.max_simulation_time


jg.runtime.ode_dahlquist_lambda1 = 1j
jg.runtime.ode_dahlquist_lambda2 = 0.1j
jg.runtime.ode_dahlquist_lambda3 = 0.01j

#####################################################
#####################################################
#####################################################

ref_ts_size = jg.runtime.max_simulation_time
ref_ts_method = f"EXP(direct)"


ts_order = 2
o=f"order={ts_order}"
ts_methods = [
            #f"ERK(ADDT(l1,l2,l3),{o})",
            #f"SS(SS(ERK(l3,{o}),ERK(l2,{o}),{o}),ERK(l1,{o}),{o})",
            #f"ETDRK(EXP(l1),ADDT(l2,l3),{o})",

            # The SL part in l2 simply uses exponential integration instead of SL
            f"SETTLS(l1,l2({o}),l3,{o})",
            f"SLETDRK(EXP(l1),l2({o}),l3,{o},version=default)",
    ]


timestep_size_min = 2**-12
timestep_sizes = [timestep_size_min*(2.0**i) for i in range(0, 6)]

#####################################################
#####################################################
#####################################################

jgbase = copy.deepcopy(jg)
jg = None

jgref = copy.deepcopy(jgbase)

#
# Reference solution
#
jgbase.runtime.rexi_method = None
jgbase.runtime.timestepping_method = ref_ts_method
jgbase.runtime.timestep_size = ref_ts_size

jgbase.reference_job = True
jgbase.gen_jobscript_directory()
jgbase.reference_job = False

# Use this one as the reference solution!
jgbase.reference_job_unique_id = jgbase.job_unique_id



#
# Create job scripts
#
for tsm in ts_methods:
    for timestep_size in timestep_sizes:

        jg = copy.deepcopy(jgbase)

        jg.runtime.timestepping_method = tsm
        jg.runtime.timestep_size = timestep_size

        if 'REXI' in jg.runtime.timestepping_method:

            if 1:
                # CI REXI via file input
                jg.runtime.rexi_method = 'file'

                #
                # !!! WARNING !!!
                # A timestep size of 32 leads to issues with the REXI method on Gitlab CI.
                # Probably, the solver accuracy is degraded due to the accummulation of the REXI terms
                #
                jg.runtime.timestep_size *= 2

                cirexi = CIREXI()
                N = 16
                R = 1
                coeffs_phi0 = cirexi.setup("phi0", N=N, R=R).toFloat()
                coeffs_phi1 = cirexi.setup("phi1", N=N, R=R).toFloat()
                coeffs_phi2 = cirexi.setup("phi2", N=N, R=R).toFloat()
                coeffs_phi3 = cirexi.setup("phi3", N=N, R=R).toFloat()

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


