#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

echo
echo "SPHERICAL HARMONICS OUTPUT"
echo

SCONS="scons --program=programs/vis_sphericalHarmonics --gui=disable --sphere2d-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit


mule.benchmark.cleanup_all || exit 1
