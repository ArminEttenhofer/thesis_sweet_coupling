#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

SCONS="scons --program=burgers --parareal=serial --cart2d-spectral-space=enable --mode=debug"
echo $SCONS
$SCONS || exit

mule.benchmark.cleanup_all || exit 1
