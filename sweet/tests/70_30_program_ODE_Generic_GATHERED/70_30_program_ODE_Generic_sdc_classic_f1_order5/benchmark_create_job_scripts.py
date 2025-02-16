#! /usr/bin/env python3

import sys
import copy
from mule.parHelper import *

from mule.JobMule import *
jg = JobGeneration()

from mule.sdc import getSDCSetup



jg.unique_id_filter = [
        "compile",
        "parallelization",
        "runtime.benchmark",
        "runtime.ode_generic",
    ]

jg.compile.program = "programs/ODE_Generic"
jg.compile.mode = "release"


jg.runtime.verbosity = 2
jg.runtime.compute_errors = 0
jg.runtime.benchmark_name = "initone"
jg.runtime.ode = "dahlquist"

jg.runtime.max_simulation_time = 10

jg.runtime.output_timestep_size = jg.runtime.max_simulation_time


ts_order = 5

nNodes = 4

paramsSDC_ = []

paramsSDC_ += [getSDCSetup(
    nNodes=nNodes,
    nIter=ts_order-1,
    nodeType="RADAU-RIGHT",
    postSweep="LASTNODE"
)]

if ts_order >= 2:
    paramsSDC_ += [getSDCSetup(
        nNodes=nNodes,
        nIter=ts_order-2,
        nodeType="RADAU-RIGHT",
        postSweep="QUADRATURE"
    )]

if ts_order >= 2:
    paramsSDC_ += [getSDCSetup(
        nNodes=nNodes,
        nIter=ts_order-2,
        nodeType="RADAU-LEFT",
        postSweep="QUADRATURE"
    )]


#####################################################
#####################################################
#####################################################

ref_ts_size = jg.runtime.max_simulation_time
ref_ts_method = f"EXP(direct)"


ts_methods = [
        f"SDCClassic(e=f)",
        f"SDCClassic(i=f)",
        f"SDCClassic(ti=ERK(f,o=1),tit=f)",
    ]


jg.runtime.ode_dahlquist_lambda1 = 0.0j
jg.runtime.ode_dahlquist_lambda2 = 0.0j
jg.runtime.ode_dahlquist_lambda3 = 0.0
jg.runtime.ode_dahlquist_mu = 0
jg.runtime.ode_dahlquist_mu = 100
jg.runtime.ode_dahlquist_phi = 20.0j
#jg.runtime.ode_dahlquist_phi = 0


timestep_size_min = 2**-12
timestep_sizes = [timestep_size_min*(2.0**i) for i in range(4, 10)]
#timestep_sizes = [timestep_size_min*(2.0**i) for i in range(0, 8)]

#####################################################
#####################################################
#####################################################


jgbase = copy.deepcopy(jg)
jg = None

jgref = copy.deepcopy(jgbase)

#
# Reference solution
#
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

        if "SDC" in tsm:
            for paramsSDC in paramsSDC_:

                jg = copy.deepcopy(jgbase)

                jg.runtime.paramsSDC = paramsSDC

                jg.runtime.timestepping_method = tsm
                jg.runtime.timestep_size = timestep_size

                setupParallelization(jg)
                jg.gen_jobscript_directory()
                
        else:

                jg = copy.deepcopy(jgbase)

                jg.runtime.timestepping_method = tsm
                jg.runtime.timestep_size = timestep_size

                setupParallelization(jg)
                jg.gen_jobscript_directory()
