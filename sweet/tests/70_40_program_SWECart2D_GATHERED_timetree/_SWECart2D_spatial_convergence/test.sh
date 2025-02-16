#! /bin/bash

#######cd "$(dirname $0)"
#######
#######mule.benchmark.cleanup_all || exit 1
#######
#######./benchmark_create_jobs.py || exit 1
#######
#######mule.benchmark.jobs_run_directly || exit 1
#######
#######./postprocessing_convergence_test.py || exit 1
#######
#######mule.benchmark.cleanup_all || exit 1


cd "$(dirname $0)"


TIMESTEPPING_GROUP="ln2space"

COMMON="../SWECart2D_timestepper_convergence_common_no_test/"


mule.benchmark.cleanup_all || exit 1

$COMMON/benchmark_create_job_scripts.py $TIMESTEPPING_GROUP || exit 1

mule.benchmark.jobs_run_directly || exit 1

$COMMON/postprocessing_pickle.py || exit 1

$COMMON/postprocessing_convergence_test.py || exit 1

mule.benchmark.cleanup_all || exit 1
