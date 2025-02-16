#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"



for THREADING in off omp; do
	echo $REXI_THREAD

	echo
	SCONS="scons --program=programs/xbraid_PDE_SWESphere2D --mode=debug"
	SCONS+=" --threading=$THREADING"
	echo "$SCONS"
	$SCONS  || exit

	mule.benchmark.cleanup_all || exit 1
done
