#! /usr/bin/env python3

#
# Delete this file if you find me. Has just testing purpose!
#

import sys
import os
import copy
import numpy as np


d = os.path.dirname(os.path.realpath(__file__))
sys.path.append(d+"/..")

from Functions import *
from cirexi.CIREXI import *

sys.path.pop()

rexi_method = "cirexi"
function_name = "phi0"


# CI-REXI: Number of quadrature poles
N = 16

# CI-REXI parameters
ci_rexi_params = {
        #'lambda_include_imag': 2,   # CI-REXI: Value on imaginary axis to be included
        #'lambda_max_real': 2        # CI-REXI: Maximum value of quadrature pole
        'R': np.pi, # radius or y extension
        'Rx': 0.3 # x extension
    }


# number of concatenations
Nconcat = int(105/ci_rexi_params['R'])

# Testing: number of samples
num_test_samples = 80

# Testing: Range (start, end)
test_range = [0, 100.0]

# Error to test for
error_eps = 1e-8

# efloat_mode
efloat_mode = "float"


def rexiConcatenate(coeffs1, coeffs2):
    """
    Compute new beta coefficients
    """
    M = len(coeffs1.alphas)
    N = len(coeffs2.alphas)

    betasp1 = [None for _ in range(M)]
    for m in range(M):
        betasp1[m] = 0
        for n in range(N):
            try:
                betasp1[m] += coeffs2.betas[n]*coeffs1.betas[m]/(coeffs1.alphas[m] - coeffs2.alphas[n])
            except:
                print("ERROR: ", coeffs1.alphas[m], coeffs2.alphas[n])
                raise Exception("ERROR")

    betasp2 = [None for _ in range(N)]
    for n in range(N):
        betasp2[n] = 0
        for m in range(M):
            try:
                betasp2[n] += -coeffs2.betas[n]*coeffs1.betas[m]/(coeffs1.alphas[m] - coeffs2.alphas[n])
            except:
                print("ERROR: ", coeffs1.alphas[m], coeffs2.alphas[n])
                raise Exception("ERROR")

    """
    Merge coefficients
    """
    coeffs = copy.deepcopy(coeffs1)

    coeffs.alphas += coeffs2.alphas
    coeffs.betas = betasp1+betasp2

    return coeffs


for m in [1]:

    if m == 0:
        print("*"*80)
        print("Classical")
        print("*"*80)

        dt_factor = 1

        Nconcat
        cirexi = CIREXI(efloat_mode=efloat_mode)
        coeffs = cirexi.setup(
            function_name = function_name,
            N = N,
            **ci_rexi_params,
        )
        coeffs = coeffs.toFloat()

        print("")
        print("alphas")
        for a in coeffs.alphas:
            print(f"    {a}")

        print("")
        print("betas")
        for b in coeffs.betas:
            print(f"    {b}")


    elif m == 1:
        print("*"*80)
        print("REXI 3000")
        print("*"*80)

        dt_factor = Nconcat

        coeffs = None

        for i in range(Nconcat):
            cirexi = CIREXI(efloat_mode=efloat_mode)
            coeffs_ = cirexi.setup(
                function_name = function_name,
                N = N,
                **ci_rexi_params,
                pole_shift_value = i/Nconcat,
            )
            coeffs_ = coeffs_.toFloat()

            if i == 0:
                coeffs = coeffs_
                continue

            coeffs = rexiConcatenate(coeffs, coeffs_)


        if 1:

            print("")
            print("alphas")
            for a in coeffs.alphas:
                print(f"    {a}")

            print("")
            print("betas")
            for b in coeffs.betas:
                print(f"    {b}")


unique_id_string = cirexi.getUniqueId()


function = Functions(
    function_name = function_name,
    efloat_mode = "float"
)


print("")
print("COEFF_NR\tALPHAS\t\t\tBETAS")
for i in range(len(coeffs.alphas)):
    print(str(i)+"\t"+str(coeffs.alphas[i])+"\t"+str(coeffs.betas[i]))

max_error = 0
verbosity = 0
for x in np.linspace(test_range[0], test_range[1], num_test_samples):
    lam = 1j*x

    y = function.eval(lam)
    yn = coeffs.eval(lam/dt_factor)

    err = np.abs(y-yn)

    print(x, err)

    if verbosity > 0:
        #if True:
        if False:
            print("x="+str(lam)+"\t\terror="+str(err))
        else:
            print("Lambda: "+str(lam))
            print(" +  exact: "+str(y))
            print(" + approx: "+str(yn))
            print(" + Error: "+str(err))
            print("")

    max_error = max(max_error, err)


if verbosity == 0:
    print(" + test_range: ["+str(test_range[0])+", "+str(test_range[1])+"]")
    print(" + Error: "+str(max_error))

if max_error > error_eps:
    raise Exception("Error threshold "+str(error_eps)+" exceeded")

