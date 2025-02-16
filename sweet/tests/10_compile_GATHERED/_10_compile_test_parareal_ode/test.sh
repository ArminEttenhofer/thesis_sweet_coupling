#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

## --cart2d-spectral-space=disable->enable
echo
echo "PARAREAL ODE"
SCONS="scons --program=programs/parareal_ode --parareal=serial --gui=disable --cart2d-spectral-space=enable --mode=debug "
echo "$SCONS"
$SCONS || exit

mule.benchmark.cleanup_all || exit 1
