#! /bin/bash

cd "$(dirname $0)"

for ORDER in 1 2 3 4; do
	mule.benchmark.cleanup_all || exit 1

	./benchmark_create_job_scripts.py $ORDER || exit 1

	mule.benchmark.jobs_run_directly || exit 1

	mule.postprocessing.pickle.alljobs.generic_data_norms || exit 1

	./postprocessing_convergence_test.py $ORDER || exit 1

	mule.benchmark.cleanup_all || exit 1
done
