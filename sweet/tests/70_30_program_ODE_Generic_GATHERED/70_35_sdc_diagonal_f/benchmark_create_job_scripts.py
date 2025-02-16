#! /usr/bin/env python3
import numpy as np
from mule.sdc import getSDCSetup
from mule.sdc.testing import setupDahlquistBenchmark, setupSDCConvergenceRuns

lamI = 0.0j
lamE = 0.0j
mu = 20.0
phi = 4*np.pi*1j
tEnd = np.pi + np.pi/2

nStepMin = 32        # minimum number of time-steps per period
nRefinement = 2      # number of time-step refinements
dtSizes = [2*np.pi/(nStepMin * 2**i) for i in range(nRefinement+1)]
print("Time steps :", dtSizes)

# Setup Dahlquist benchmark and reference simulation
jg, terms = setupDahlquistBenchmark(tEnd, lamI, lamE, mu, phi, fExplicit=False)

# Setup SDC runs
baseParams = dict(
    nNodes=None,
    nIter=None,
    nodeType=None,
    qDeltaImplicit=["DNODES-3", "DNODES-4"],
    qDeltaExplicit="PIC",
    preSweep="QDELTA",
    qDeltaInitial="TRAPAR",
    postSweep="LASTNODE",
    nodeDistr=None,
)
for nIter in [0, 1, 2, 3]:
    for nNodes in [1, 2, 3, 4]:
        for nodeDistr in ['LEGENDRE', 'EQUID']:
            for nodeType in ['LOBATTO', 'RADAU-RIGHT']:
                if nodeType == 'LOBATTO' and nNodes == 1:
                    continue
                # Handle edge cases ---------------------------------------------------------------------
                expectedOrder = None
                if nNodes == 2 and nodeDistr == "LEGENDRE" and nodeType == "RADAU-RIGHT" and nIter == 1:
                    expectedOrder = 3
                if nNodes == 3:
                    if nIter == 1:
                        if nodeType == "LOBATTO":
                            expectedOrder = 4
                        if nodeType == "RADAU-RIGHT":
                            if nodeDistr == "EQUID":
                                expectedOrder = 3
                            if nodeDistr == "LEGENDRE":
                                expectedOrder = 5
                    if nIter == 2 and nodeType == "RADAU-RIGHT" and nodeDistr == "LEGENDRE":
                        expectedOrder = 5
                if nNodes == 4:
                    if nodeType == "LOBATTO":
                        if nIter == 1:
                            if nodeDistr == "EQUID":
                                expectedOrder = 4
                            if nodeDistr == "LEGENDRE":
                                expectedOrder = 6
                        if nIter in [2, 3] and nodeDistr == "LEGENDRE":
                            expectedOrder = 6
                    if nodeType == "RADAU-RIGHT":
                        if nIter == 1: 
                            if nodeDistr == "EQUID":
                                expectedOrder = 4
                            if nodeDistr == "LEGENDRE":
                                expectedOrder = 7
                        if nIter in [2, 3] and nodeDistr == "LEGENDRE":
                            expectedOrder = 7
                # ---------------------------------------------------------------------------------------
                params = baseParams.copy()
                params.update(nIter=nIter, nNodes=nNodes, nodeDistr=nodeDistr, nodeType=nodeType)
                params = getSDCSetup(**params)
                setupSDCConvergenceRuns(
                    jg, terms, dtSizes, params, formulation="FP", expectedOrder=expectedOrder)
