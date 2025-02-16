import itertools
import copy
import numpy as np

from mule.JobMule import JobGeneration
from mule.parHelper import setupParallelization
from mule.etdsdc import getETDSDCSetup
from mule.SWEETFileDict import SWEETFileDict

from mule.postprocessing.JobsData import JobsData
from mule.postprocessing.JobsDataConsolidate import JobsDataConsolidate, JobsData_GroupsPlottingScattered, JobsData_GroupsCleanupPostprocessed

# Matplotlib has to be imported after mule.postprocessing
import matplotlib.pyplot as plt

TEST_CONFIGS = {
    "nNodes": [1, 2, 3], # , 4, 5, 7, 9],
    "nodeDistr": ["CHEBY-1", "LEGENDRE", "EQUID"],
}

# get all permutations
def setupTestParams(test_configs=TEST_CONFIGS):
    # NOTE: sort in alphabetic. order to group for plotting later
    if "nodeDistr" in test_configs:
        test_configs["nodeDistr"].sort()
    if "nNodes" in test_configs:
        test_configs["nNodes"].sort()
    if "nIter" in test_configs:
        test_configs["nIter"].sort()

    test_params = []
    keys = list(test_configs.keys())
    for conf in itertools.product(*test_configs.values()):
        test_params.append({k: v for k, v in zip(keys, conf)})

    return test_params


def setupDahlquistBenchmark(tEnd, lamI=0.0, lamE=0.0, mu=0.0, phi=0.0, fExplicit=True):
    jg = JobGeneration()
    jg.unique_id_filter = [
            "compile",
            "parallelization",
            "runtime.benchmark",
            "runtime.ode_generic",
        ]

    jg.compile.program = "programs/ODE_Generic"
    jg.compile.mode = "debug"

    jg.runtime.verbosity = 2
    jg.runtime.compute_errors = 0
    jg.runtime.benchmark_name = "initone"
    jg.runtime.ode = "dahlquist"

    jg.runtime.ode_dahlquist_lambda1 = lamI
    jg.runtime.ode_dahlquist_lambda2 = lamE
    jg.runtime.ode_dahlquist_lambda3 = 0.0
    jg.runtime.ode_dahlquist_mu = mu
    jg.runtime.ode_dahlquist_phi = phi

    jg.runtime.max_simulation_time = tEnd
    jg.runtime.output_timestep_size = tEnd

    ref_ts_size = tEnd
    ref_ts_method = "EXP(direct)"  # exponential integration for reference solution

    # Reference solution
    jg.runtime.timestepping_method = ref_ts_method
    jg.runtime.timestep_size = ref_ts_size

    jg.reference_job = True
    jg.gen_jobscript_directory()
    jg.reference_job = False

    # Use this one as the reference solution!
    jg.reference_job_unique_id = jg.job_unique_id

    terms = []

    iTerm = ""
    if lamI != 0.0:
        iTerm = "l1"
    elif mu != 0 and not fExplicit:
        iTerm = "f"

    eTerm = ""
    if lamE != 0 and mu != 0:
        eTerm = "ADDT(l2,f)"
    elif lamE != 0:
        eTerm = "l2"
    elif mu != 0 and fExplicit:
        eTerm = "f"
    
    if iTerm: 
        terms.append(iTerm)
    if eTerm: 
        terms.append(eTerm)
    terms = ",".join(terms)

    return jg, terms


def setupSDCConvergenceRuns(
        jgBase:JobGeneration, terms, dtSizes, paramsSDC,
        expectedOrder=None):
    
    jgGroup = copy.deepcopy(jgBase)
    paramsSDC = paramsSDC if isinstance(paramsSDC, SWEETFileDict) else getETDSDCSetup(**paramsSDC)
    order = min(paramsSDC["nIter"] + 1, len(paramsSDC["tau"])) # M+1, N
    method = f"ETDSDC({terms})"
    jgGroup.runtime.timestepping_method = method
    jgGroup.runtime.paramsSDC = paramsSDC
    jgGroup.runtime.init_phase = True
    jgGroup.runtime.idString = method+'_'+paramsSDC['idString']
    jgGroup.runtime.sdcOrder = order if expectedOrder is None else expectedOrder
    jgGroup.runtime.init_phase = False
    
    for dt in dtSizes:
        jgRun = copy.deepcopy(jgGroup)
        jgRun.runtime.timestep_size = dt
        setupParallelization(jgRun)
        jgRun.gen_jobscript_directory()
            

def checkConvergence(test_configs=TEST_CONFIGS, groupby="None", generatePlot=True, checkOrder=False, plotLim=None):

    j = JobsData(verbosity=0)
    c = JobsDataConsolidate(j)

    # Group together SDC runs with same idString
    groups = ['runtime.idString']
    print("Groups:")
    job_groups = c.create_groups(groups)

    # Cleanup postprocessed data
    tag_cleanup_info = [
        {
        "ref_file_starts_with": "output_", 
        "tag_src": "res_norm_linf", 
        "tag_dst": "generic_data_norms.res_norm_linf"
        }
    ]
    JobsData_GroupsCleanupPostprocessed(job_groups, tag_cleanup_info, pickle_file_default_prefix="generic_data_norms_")

    # Use plotting function to create (x,y) data
    d = JobsData_GroupsPlottingScattered(
        job_groups,
        'runtime.timestep_size',
        'generic_data_norms.res_norm_linf',
        meta_attribute_name = 'runtime.sdcOrder',
    )

    for gName, gData in d.get_data_float().items():
        print("*"*80)
        print("Group: " + gName)
        
        # NOTE: while same ETDSDC is tested - shorten gName
        gName = gName[gName.rfind(")")+2:]

        prev_err = None
        conv = '-'

        dtSizes = gData['x_values']
        errors = gData['y_values']
        
        if generatePlot:
            plt.loglog(dtSizes, errors, 'o-', label=gName)
            plt.legend()

        # make sure of descending order of dt
        for dt, err in sorted(zip(dtSizes, errors), key=lambda x: -x[0]):
            if prev_err == None:
                conv = '-----'
            else:
                conv = f"{np.log2(prev_err/err):.2f}"

            expect = int(gData['meta_values'][0]) if conv != '-----' else '--'
            print(f"\t{dt:.5e}\t=>\t{err:.16e}\tconvergence: {conv}\t | expected: {expect}")
            prev_err = err

            if checkOrder and conv != "[error=0]":
                if np.isnan(conv):
                    print("ERROR: got nan, but not for a kebab !")
                    raise Exception(f"Got nan for {gName}")
                conv = max(conv, 0)
                if abs(conv-expect) > 0.3:
                    print("ERROR: not the expected order !")
                    raise Exception(f"Convergence mismatch for {gName}")

        # group all nodeDistr in 1 plot, since "LEGENDRE" is last
        if generatePlot and groupby=="nodeDistr" and (test_configs["nodeDistr"][-1] in gName):
            plt.grid()
            plt.ylim(plotLim)
            plt.xlabel(r'$\Delta{t}$')
            plt.ylabel(r'$L_\infty$ error')
            plt.title("ETDSDC: Dahlquist")
            figname = "conv_ETDSDC_" + gName + ".pdf"
            plt.savefig(figname, bbox_inches='tight')
            plt.clf()
        
    if generatePlot and groupby=="None":
        plt.grid()
        plt.ylim(plotLim)
        plt.xlabel(r'$\Delta{t}$')
        plt.ylabel(r'$L_\infty$ error')
        plt.title("ETDSDC: Dahlquist")
        figname = "conv_ETDSDC.pdf"
        plt.savefig(figname, bbox_inches='tight')