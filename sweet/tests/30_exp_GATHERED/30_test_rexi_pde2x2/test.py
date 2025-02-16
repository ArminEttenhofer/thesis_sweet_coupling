#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from mule.JobMule import *
from mule.utils import exec_program


#
# REXI specific
#
from mule.rexi.REXICoefficients import *
from mule.rexi.trexi.TREXI import *
from mule.rexi.cirexi.CIREXI import *
from mule.rexi.brexi.BREXI import *


#
# Cleanup potentially existing benchmarks
#
exec_program('mule.benchmark.cleanup_all', catch_output=False)

from itertools import product


#
# Generate new job generation instance
#
jg = JobGeneration()


jg.compile.program = "tests/exp_rexiPDE2x2"

#
# Use file-based REXI method
#
jg.runtime.rexi_method = "file"
jg.runtime.verbosity = 20


rexi_file_methods = [
        "trexi",
        "cirexi",
        "brexi"
    ]


#function_name_list = ["phi0", "phi1", "phi2", "phi3", "phi4", "ups1", "ups2", "ups3", "psi1", "psi2", "psi3"]
function_name_list = ["phi0", "phi1", "phi2", "phi3", "phi4", "ups1", "ups2", "ups3"]


efloat_mode = "float"
#efloat_mode = "mpfloat"



for function_name in function_name_list:

    for rexi_file_method in rexi_file_methods:

        if rexi_file_method == "trexi":

            if function_name not in ["phi0", "phi1", "phi2"]:
                continue

            trexi = TREXI(efloat_mode = efloat_mode)

            if function_name in ["phi1", "phi2"]:
                M_list = [128, 256]
                h_list = [0.1, 0.2, 0.5]
            else:
                M_list = [64, 128, 256]
                h_list = [0.1, 0.2, 0.5]

            for (M, h) in product(M_list, h_list):
                coeffs = trexi.setup(M=M, h=h).toFloat()

                if function_name == "phi1":
                    coeffs = coeffs.getNextOrderPhi()

                elif function_name == "phi2":
                    coeffs = coeffs.getNextOrderPhi()
                    coeffs = coeffs.getNextOrderPhi()

                coeffs.normalize_steady_state()

                jg.runtime.rexi_files_coefficients = [coeffs]
                jg.gen_jobscript_directory()

                if True:
                    # Validate with C-implementation of T-REXI method
                    jg.runtime.rexi_method = "terry"
                    jg.runtime.rexi_terry_m = M
                    jg.runtime.rexi_terry_h = h
                    jg.gen_jobscript_directory()

                    # Back to original version
                    jg.runtime.rexi_method = "file"


        elif rexi_file_method == "cirexi":

            cirexi = CIREXI(efloat_mode = efloat_mode)

            # CI-REXI: Number of quadrature poles
            N_list = [256, 512]

            # CI-REXI: Value on imaginary axis to be included
            lambda_include_imag_list = [15, 20]

            # CI-REXI: Maximum value of quadrature pole
            lambda_max_real_list = [5, 10]

            for (N, lambda_include_imag, lambda_max_real) in product(N_list, lambda_include_imag_list, lambda_max_real_list):
                coeffs = cirexi.setup(function_name=function_name, N=N, lambda_include_imag=lambda_include_imag, lambda_max_real=lambda_max_real).toFloat()
                jg.runtime.rexi_files_coefficients = [coeffs]
                jg.gen_jobscript_directory()

                if True:
                    # Validate with C-implementation of CI-REXI method
                    jg.runtime.rexi_method = "ci"
                    jg.runtime.rexi_ci_n = N
                    jg.runtime.rexi_ci_max_real = lambda_max_real
                    jg.runtime.rexi_ci_max_imag = lambda_include_imag
                    jg.gen_jobscript_directory()

                    # Back to original version
                    jg.runtime.rexi_method = "file"


        elif rexi_file_method == "brexi":

            if function_name not in ["phi0", "phi1", "phi2"]:
                continue

            brexi = BREXI(efloat_mode = efloat_mode)

            N_list = [8, 10]
            N_list = [12]

            quadrature_method_list = ["gauss_legendre", "gauss_chebyshev_u"]
            quadrature_method_list = ["gauss_chebyshev_u"]
            #quadrature_method_list = ["gauss_legendre"]

            for quadrature_method in quadrature_method_list:

                #
                # We use quite large time steps and need to fine tune the number of quadrature poles
                #

                if quadrature_method == "gauss_legendre":
                    if function_name == "phi0":
                        N_list = [6]
                    elif function_name == "phi1":
                        N_list = [12]
                    elif function_name == "phi2":
                        N_list = [12]
                    else:
                        raise Exception("TODO")

                elif quadrature_method == "gauss_chebyshev_u":
                    if function_name == "phi0":
                        N_list = [6]
                    elif function_name == "phi1":
                        N_list = [16]
                    elif function_name == "phi2":
                        N_list = [16]
                    else:
                        raise Exception("TODO")

                else:
                    raise Exception("TODO")

                for N in N_list:

                    coeffs = brexi.setup(N=N, quadrature_method=quadrature_method).toFloat()

                    if function_name == "phi1":
                        coeffs = coeffs.getNextOrderPhi()

                    elif function_name == "phi2":
                        coeffs = coeffs.getNextOrderPhi()
                        coeffs = coeffs.getNextOrderPhi()

                    coeffs = coeffs.toFloat()

                    jg.runtime.rexi_files_coefficients = [coeffs]
                    jg.gen_jobscript_directory()

        else:
            raise Exception("Unknown method "+rexi_file_methods)



exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
    sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)
