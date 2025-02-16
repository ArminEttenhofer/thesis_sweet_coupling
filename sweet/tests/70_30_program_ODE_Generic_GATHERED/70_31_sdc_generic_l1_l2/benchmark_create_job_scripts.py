#! /usr/bin/env python3
import numpy as np
from mule.sdc import getSDCSetup
from mule.sdc.testing import setupDahlquistBenchmark, setupSDCConvergenceRuns, TEST_PARAMS

lamI = 1.0j
lamE = 0.1j
mu = 0.0
tEnd = np.pi + np.pi/2

def genTimeSteps(nStepMin=64, nRefinement=2):
    return [2*np.pi/(nStepMin*2**i) for i in range(nRefinement+1)]

# Setup Dahlquist benchmark and reference simulation
jg, terms = setupDahlquistBenchmark(tEnd, lamI, lamE, mu)

for params in TEST_PARAMS:
    params = getSDCSetup(**params)
    dtSizes = genTimeSteps()
    setupSDCConvergenceRuns(jg, terms, dtSizes, params, formulation="FP")
