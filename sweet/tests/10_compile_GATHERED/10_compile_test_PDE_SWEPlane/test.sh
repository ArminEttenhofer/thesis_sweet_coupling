#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

SCONS="scons --program=programs/PDE_SWECart2D --gui=disable --cart2d-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit

####echo
####echo "SWE REXI"
####SCONS="scons --program=programs/PDE_SWECart2D --gui=disable --cart2d-spectral-space=disable --libfft=enable --mode=debug"
####echo "$SCONS"
####$SCONS  || exit

###echo
###echo "SWE REXI"
###SCONS="scons --program=programs/PDE_SWECart2D --gui=disable --cart2d-spectral-space=disable --libfft=enable --mode=debug"
###echo "$SCONS"
###$SCONS  || exit


SCONS="scons --program=programs/PDE_SWECart2D --gui=disable --cart2d-spectral-space=enable --cart2d-spectral-dealiasing=disable --mode=debug"
echo "$SCONS"
$SCONS  || exit

if [ "$SWEET_MPICXX" != "" ]; then
	SCONS="scons --program=programs/PDE_SWECart2D --sweet-mpi=enable --rexi-thread-parallel-sum=enable --threading=off"
	echo "$SCONS"
	$SCONS  || exit
fi

mule.benchmark.cleanup_all || exit 1
