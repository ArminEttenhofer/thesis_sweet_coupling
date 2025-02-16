#! /usr/bin/env python3
import numpy as np
from mule.sdc import getSDCSetup
from mule.sdc.testing import setupDahlquistBenchmark, setupSDCConvergenceRuns

lamI = 1.0j
lamE = 0.0j
mu = 0.0
tEnd = np.pi + np.pi/2

nStepMin = 64        # minimum number of time-steps per period
nRefinement = 2      # number of time-step refinements
dtSizes = [2*np.pi/(nStepMin * 2**i) for i in range(nRefinement+1)]
print("Time steps :", dtSizes)

# Setup Dahlquist benchmark and reference simulation
jg, terms = setupDahlquistBenchmark(tEnd, lamI, lamE, mu)

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
                params = baseParams.copy()
                params.update(nIter=nIter, nNodes=nNodes, nodeDistr=nodeDistr, nodeType=nodeType)
                params = getSDCSetup(**params)
                setupSDCConvergenceRuns(jg, terms, dtSizes, params, formulation="FP")
