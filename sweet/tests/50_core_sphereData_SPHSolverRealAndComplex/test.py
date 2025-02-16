#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from mule.JobMule import *
from itertools import product
from mule.utils import exec_program

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = JobGeneration()
jg.compile.program="tests/core_sphere2dDataSPHSolverRealAndComplex"

jg.compile.cart2d_spectral_space="disable"
jg.compile.sphere2d_spectral_space="enable"
jg.compile.mode = "release"

# Request dedicated compile script
jg.compilecommand_in_jobscript = False


jg.runtime.sphere2d_radius = 1
jg.runtime.sphere2d_rotating_coriolis_omega = 1

unique_id_filter = []
unique_id_filter.append('compile')


jg.unique_id_filter = unique_id_filter


#params_runtime_mode_res = [64, 128, 256, 512, 1024, 2048]
params_runtime_mode_res = [64, 128, 256, 512, 1024]

params_runtime_r = [1, 1e3, 1e6]
params_runtime_f = [1, 1e-3, 1e-6]

jg.runtime.verbosity = 5

for (
    	jg.runtime.space_res_spectral,
    	jg.runtime.sphere2d_radius,
    	jg.runtime.sphere2d_rotating_coriolis_omega,
    ) in product(
    	params_runtime_mode_res,
    	params_runtime_r,
    	params_runtime_f,
    ):
    jg.gen_jobscript_directory()


jg.write_compilecommands(compilecommands_filename="./compile_platform.sh")

(output, exitcode) = exec_program('./compile_platform.sh', catch_output=True)
if exitcode != 0:
    print("Output:")
    print(output)
    sys.exit(exitcode)

(output, exitcode) = exec_program('mule.benchmark.jobs_run_directly', catch_output=True)
if exitcode != 0:
    print("Output:")
    print(output)
    sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)

