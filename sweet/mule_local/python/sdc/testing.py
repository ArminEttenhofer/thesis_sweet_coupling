import itertools
import copy
import numpy as np
import re

from mule.JobMule import JobGeneration
from mule.parHelper import setupParallelization
from mule.sdc import getSDCSetup
from mule.SWEETFileDict import SWEETFileDict

from mule.postprocessing.JobsData import JobsData
from mule.postprocessing.JobsDataConsolidate import JobsDataConsolidate, JobsData_GroupsPlottingScattered, JobsData_GroupsCleanupPostprocessed

# Matplotlib has to be imported after mule.postprocessing
import matplotlib.pyplot as plt

TEST_CONFIGS = {
    "nIter": [0, 1, 2],
    "nNodes": [1, 2, 3, 4],
    "nodeType": ["GAUSS", "LOBATTO", "RADAU-RIGHT", "RADAU-LEFT"],
    "nodeDistr": ["LEGENDRE", "EQUID"],
    "qDeltaImplicit": ["BE"],
    "qDeltaExplicit": ["FE"],
    "qDeltaInitial": ["BE"],
    "preSweep": ["QDELTA", "COPY", "ZEROS"],
    "postSweep": ["LASTNODE", "QUADRATURE", "INTERPOLATION"],
}

TEST_PARAMS = []
keys = list(TEST_CONFIGS.keys())
for conf in itertools.product(*TEST_CONFIGS.values()):
    if conf[2] in ["GAUSS", "RADAU-LEFT"] and conf[-1] == "LASTNODE":
        continue
    if conf[2] in ["LOBATTO", "RADAU-LEFT"] and conf[1] == 1:
        continue
    TEST_PARAMS.append({k: v for k, v in zip(keys, conf)})


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

    #
    # Reference solution
    #
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
        iTerm = "i=l1"
    elif mu != 0 and not fExplicit:
        iTerm = "i=f"

    eTerm = ""
    if lamE != 0 and mu != 0:
        eTerm = "e=ADDT(l2,f)"
    elif lamE != 0:
        eTerm = "e=l2"
    elif mu != 0 and fExplicit:
        eTerm = "e=f"
    
    if iTerm: 
        terms.append(iTerm)
    if eTerm: 
        terms.append(eTerm)
    terms = ",".join(terms)

    return jg, terms


def setupSDCConvergenceRuns(
        jgBase:JobGeneration, terms, dtSizes, paramsSDC, formulation="FP",
        expectedOrder=None):
    
    jgGroup = copy.deepcopy(jgBase)
    paramsSDC = paramsSDC if isinstance(paramsSDC, SWEETFileDict) else getSDCSetup(**paramsSDC)
    
    if 'i=' in terms and 'e=' in terms:
        order = min(paramsSDC['orderI'], paramsSDC['orderE'])
    elif 'i=' in terms:
        order = paramsSDC['orderI']
    elif 'e=' in terms:
        order = paramsSDC['orderE']
    else:
        raise ValueError(f'problem here with terms={terms}')

    method = f"SDC{formulation}({terms})"
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
            

def checkConvergence(generatePlot=True, checkOrder=False, plotLim=None):

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
        print("Group: "+gName)

        prev_err = None
        conv = '-'

        dtSizes = gData['x_values']
        errors = gData['y_values']
        nIters = re.search(r'(?<=_K)\d+', gName).group(0)
        
        if generatePlot:
            plt.loglog(dtSizes, errors, 'o-', label=f"nSweep={nIters}")
            plt.legend()

        for dt, err in zip(dtSizes, errors):

            if prev_err == None:
                conv = '[error=0]'
            else:
                conv = np.log2(err/prev_err)

            expect = int(gData['meta_values'][0]) if conv != '[error=0]' else '[error=0]'
            print(f"\t{dt}\t=>\t{err}\tconvergence: {conv}\t | expected: {expect}")
            prev_err = err

            if checkOrder and conv != "[error=0]":
                if np.isnan(conv):
                    print("ERROR: got nan, but not for a kebab !")
                    raise Exception(f"Got nan for {gName}")
                conv = max(conv, 0)
                if abs(conv-expect) > 0.3:
                    print("ERROR: not the expected order !")
                    raise Exception(f"Convergence mismatch for {gName}")

    if generatePlot:
        plt.grid()
        plt.ylim(plotLim)
        plt.xlabel(r'$\Delta{t}$')
        plt.ylabel(r'$L_\inf$ error')
        plt.savefig('convergence.pdf', bbox_inches='tight')