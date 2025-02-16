#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

pwd

for i in src/tests/*.cpp; do

	SRC=$(basename $i)
	UNIT_TEST="${SRC/.cpp/}"

	SCONS="scons --gui=disable --sphere2d-spectral-space=enable --mode=debug"
	SCONS+=" --cart2d-spectral-space=enable"
	SCONS+=" --sphere2d-spectral-space=enable"
	SCONS+=" --threading=omp"
	SCONS+=" --eigen=enable"
	SCONS+=" --parareal=serial"
	SCONS+=" --parareal-sphere2d=enable"
	SCONS+=" --program=tests/$UNIT_TEST"

	echo "$SCONS"
	$SCONS  || exit

	mule.benchmark.cleanup_all || exit 1
done
