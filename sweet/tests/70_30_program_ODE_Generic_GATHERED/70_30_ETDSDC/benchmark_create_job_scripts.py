#! /usr/bin/env python3
import copy
from mule.parHelper import *
from mule.JobMule import *
from mule.etdsdc import getETDSDCSetup
### compare with IMEXSDC ###
from mule.sdc import getSDCSetup


lam1 = 10.0j
lam2 = -0.03j
mu = 1e-3
phi = 1.0j

tEnd = 10. # np.pi + np.pi/2
timestep_sizes = [2.0**(-i) for i in range(2, 14)]


############################

jg = JobGeneration()

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

jg.runtime.ode_dahlquist_lambda1 = lam1
jg.runtime.ode_dahlquist_lambda2 = lam2
jg.runtime.ode_dahlquist_lambda3 = 0.0
jg.runtime.ode_dahlquist_mu = mu
jg.runtime.ode_dahlquist_phi = phi

jg.runtime.max_simulation_time = tEnd
jg.runtime.output_timestep_size = jg.runtime.max_simulation_time

ref_ts_method = "EXP(direct)"

ts_methods = ["ETDRK(l1,ADDT(l2,f),o=1)", \
              "ETDRK(l1,ADDT(l2,f),o=2)", "ETDRK(l1,ADDT(l2,f),o=4)", \
              "LRK(l1,ADDT(l2,f),o=1)", "LRK(l1,ADDT(l2,f),o=2)", \
              "LRK(l1,ADDT(l2,f),o=3)", "LRK(l1,ADDT(l2,f),o=4)", \
              "ETDSDC(l1,ADDT(l2,f))"] # add IMEX if needed "SDC_Classic(l1,ADDT(l2,f))"

jg_base = copy.deepcopy(jg)
jg = None

for tsm in ts_methods:
    # set up ref job
    jgref = copy.deepcopy(jg_base)
    jgref.runtime.timestepping_method = ref_ts_method
    jgref.runtime.timestep_size = timestep_sizes[0] # it's direct

    jgref.reference_job = True
    jgref.gen_jobscript_directory()
    jgref.reference_job = False

    #### for ETDSDC generate params ####
    if (tsm.find("ETDSDC") != -1):
        for order in range(1, 5):
            SDCparams= getETDSDCSetup(nNodes=max(order-2, 0), nIter=max(order-1, 0), nodeDistr="LEGENDRE")

            jgSDC = copy.deepcopy(jg_base)
            jgSDC.runtime.timestepping_method = tsm
            jgSDC.runtime.paramsSDC = SDCparams
            jgSDC.runtime.init_phase = True
            jgSDC.runtime.idString = tsm +'_'+SDCparams['idString']
            jgSDC.runtime.sdcOrder = order
            jgSDC.runtime.init_phase = False

            for timestep_size in timestep_sizes:
                jg = copy.deepcopy(jgSDC)

                jg.reference_job_unique_id = jgref.job_unique_id
                jg.runtime.timestep_size = timestep_size

                setupParallelization(jg)
                jg.gen_jobscript_directory()
    #### for IMEXSDC generate params ####
    elif (tsm.find("Classic") != -1):
        for order in range(1, 5):
            SDCparams= getSDCSetup(nNodes=max(order, 2), nIter=order-1, nodeDistr="LEGENDRE")

            jgSDC = copy.deepcopy(jg_base)
            jgSDC.runtime.timestepping_method = tsm
            jgSDC.runtime.paramsSDC = SDCparams
            jgSDC.runtime.init_phase = True
            jgSDC.runtime.idString = tsm +'_'+SDCparams['idString']
            jgSDC.runtime.sdcOrder = order
            jgSDC.runtime.init_phase = False

            for timestep_size in timestep_sizes:
                jg = copy.deepcopy(jgSDC)

                jg.reference_job_unique_id = jgref.job_unique_id
                jg.runtime.timestep_size = timestep_size

                setupParallelization(jg)
                jg.gen_jobscript_directory()

    #### for ETDRK, LRK just run ####
    else:
        for timestep_size in timestep_sizes:
            jg = copy.deepcopy(jg_base)

            jg.reference_job_unique_id = jgref.job_unique_id
            jg.runtime.timestepping_method = tsm
            jg.runtime.timestep_size = timestep_size
            jg.runtime.init_phase = True
            jg.runtime.idString = tsm
            jg.runtime.init_phase = False

            setupParallelization(jg)
            jg.gen_jobscript_directory()
