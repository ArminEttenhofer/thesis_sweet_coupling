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
from mule.rexi.Functions import *


jg.unique_id_filter = [
        "compile",
        "parallelization",
        "runtime.benchmark",
        "runtime.ode_generic",
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
jg.runtime.ode_dahlquist_lambda3 = 0.0

#####################################################
#####################################################
#####################################################

ref_ts_size = jg.runtime.max_simulation_time
ref_ts_method = f"EXP(direct)"


ts_order = 2
o=f"order={ts_order}"
ts_methods = [
            f"ERK(ADDT(l1,l2),{o})",

            f"SS(ERK(l1,{o}),IRK(l2,{o}),{o})",
            f"SS(EXP(l1),ERK(l2,{o}),{o})",
            f"SS(REXI(l1),ERK(l2,{o}),{o})",

            # ETDRK
            f"ETDRK(EXP(l1),l2,{o})",

            # ETDRK - REXI version
            f"ETDRK(REXI(l1),l2,{o})",

            # Test subcycling
            f"SS(SUB(ERK(l1,{o}),n=10),ERK(l2,{o}),{o})",

            # Test AddTendencies / NegTendencies
            f"SS(ERK(ADDT(l1,l2,NEGT(l2)),{o}),ERK(l2,{o}),{o})",

            # Test AddIntegration/NegIntegration
            f"SS(ADDI(ERK(l1,order=2),ERK(l2,order=2),NEGI(ERK(l2,order=2))),ERK(l2,order=2),order=2)",
    ]


timestep_size_min = 2**-12
timestep_sizes = [timestep_size_min*(2.0**i) for i in range(2, 8)]

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


        if 'rexi' not in jg.runtime.timestepping_method.lower():
            setupParallelization(jg)
            jg.gen_jobscript_directory()
            continue
        
        else:
            for method in ["cirexi", "brexi"]:

                if method == "cirexi":
                    #
                    # CI-REXI via file input
                    #
                    jg.runtime.rexi_method = 'file'

                    #
                    # !!! WARNING !!!
                    # A timestep size of 32 leads to issues with the REXI method on Gitlab CI.
                    # Probably, the solver accuracy is degraded due to the accummulation of the REXI terms
                    #
                    jg.runtime.timestep_size *= 2
                    jg.runtime.timestep_size *= 2

                    cirexi = CIREXI()
                    N = 32
                    R = 2
                    coeffs_phi0 = cirexi.setup("phi0", N=N, R=R).toFloat()

                    for promo in [True, False]:
                        if promo:
                            coeffs_phi0.unique_id_string += "_promo"
                            coeffs_phi1 = coeffs_phi0.getNextOrderPhi()
                            coeffs_phi2 = coeffs_phi1.getNextOrderPhi()

                            coeffs_phi0.normalize_steady_state()
                            coeffs_phi1.normalize_steady_state()
                            coeffs_phi2.normalize_steady_state()

                            jg.runtime.rexi_files_coefficients = [coeffs_phi0, coeffs_phi1, coeffs_phi2]

                            setupParallelization(jg)
                            jg.gen_jobscript_directory()


                        else:
                            coeffs_phi1 = cirexi.setup("phi1", N=N, R=R).toFloat()
                            coeffs_phi2 = cirexi.setup("phi2", N=N, R=R).toFloat()

                            coeffs_phi0.normalize_steady_state()
                            coeffs_phi1.normalize_steady_state()
                            coeffs_phi2.normalize_steady_state()

                            jg.runtime.rexi_files_coefficients = [coeffs_phi0, coeffs_phi1, coeffs_phi2]

                            setupParallelization(jg)
                            jg.gen_jobscript_directory()


                elif method == "brexi":
                    #
                    # B-REXI via file
                    #
                    jg.runtime.rexi_method = 'file'

                    brexi_phi0 = BREXI()
                    #brexi_phi0.setup(N=8, quadrature_method='gauss_legendre')
                    #brexi_phi0.setup(N=8, quadrature_method='gauss_jacobi')
                    brexi_phi0.setup(N=8, quadrature_method='gauss_chebyshev_u')

                    coeffs_phi0 = brexi_phi0.getCoeffs().toFloat()
                    coeffs_phi1 = coeffs_phi0.getNextOrderPhi()
                    coeffs_phi2 = coeffs_phi1.getNextOrderPhi()

                    jg.runtime.rexi_files_coefficients = [coeffs_phi0, coeffs_phi1, coeffs_phi2]

                    setupParallelization(jg)
                    jg.gen_jobscript_directory()

                else:
                    raise Exception("TODO")

