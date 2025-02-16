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


if len(sys.argv) > 1:
    ts_order = int(sys.argv[1])
else:
    ts_order = 2


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


o=f"order={ts_order}"
ts_methods = [
        #f"SDC-FP(e=f,l=l1)",
        #f"SDC-FP(i=f,n=l1)",
        #f"SDC-FP(i=l1,e=f)",
        #f"SDC-FP(i=f,e=l1)",
        f"SDCClassic(i=f,e=ADDT(l1,l2))",
        f"SDCClassic(i=l1,e=ADDT(f,l2))",
        f"SDCClassic(e=ADDT(l1,l2,f))",
        f"SDCClassic(e=f,ti=ERK(ADDT(l1,l2),o=1),tit=ADDT(l1,l2))",
        f"SDCClassic(ti=ERK(ADDT(l1,l2,f),o=1),tit=ADDT(l1,l2,f))",
        f"SDCClassic(e=ADDT(l1,l2),ti=ERK(f,o=1),tit=ADDT(f))",
    ]


if ts_order <= 4:
    ts_methods += [f"ERK(ADDT(l1,l2,f),{o})"]


jg.runtime.ode_dahlquist_lambda1 = 5.0j
jg.runtime.ode_dahlquist_lambda2 = 0.1j
jg.runtime.ode_dahlquist_lambda3 = 0.0
jg.runtime.ode_dahlquist_mu = 0
jg.runtime.ode_dahlquist_mu = 100
jg.runtime.ode_dahlquist_phi = 20.0j
#jg.runtime.ode_dahlquist_phi = 0


timestep_size_min = 2**-12
timestep_sizes = [timestep_size_min*(2.0**i) for i in range(0+ts_order, 10+ts_order)]

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
