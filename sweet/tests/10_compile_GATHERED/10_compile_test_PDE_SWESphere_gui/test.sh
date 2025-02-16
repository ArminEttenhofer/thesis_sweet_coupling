#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"


for REXI_THREAD in enable disable; do
	for SWEET_MPI in enable disable; do
		for THREADING in off omp; do
			echo $REXI_THREAD

			echo
			SCONS="scons --program=programs/PDE_SWESphere2D --mode=debug"
			SCONS+=" --gui=enable"
			SCONS+=" --sphere2d-spectral-space=enable"
			SCONS+=" --rexi-thread-parallel-sum=$REXI_THREAD"
			SCONS+=" --sweet-mpi=$SWEET_MPI"
			SCONS+=" --threading=$THREADING"
			echo "$SCONS"
			$SCONS  || exit

			mule.benchmark.cleanup_all || exit 1
		done
	done
done
