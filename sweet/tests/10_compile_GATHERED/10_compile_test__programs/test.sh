#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

pwd

for i in src/programs/[^_]*.cpp; do
	#for MPI in enable disable; do
	# Not all programs support MPI, hence it's disabled per default
	for MPI in disable; do

		SRC=$(basename $i)
		UNIT_TEST="${SRC/.cpp/}"

		SCONS="scons --gui=disable --sphere2d-spectral-space=enable --mode=debug"
		SCONS+=" --cart2d-spectral-space=enable"
		SCONS+=" --sphere2d-spectral-space=enable"
		SCONS+=" --threading=omp"
		SCONS+=" --eigen=enable"
		SCONS+=" --fortran-source=enable"
		SCONS+=" --sweet-mpi=${MPI}"
		#SCONS+=" --libpfasst=enable"
		#SCONS+=" --parareal=serial"
		SCONS+=" --program=programs/$UNIT_TEST"

		echo "$SCONS"
		$SCONS  || exit

	done
done

mule.benchmark.cleanup_all || exit 1
