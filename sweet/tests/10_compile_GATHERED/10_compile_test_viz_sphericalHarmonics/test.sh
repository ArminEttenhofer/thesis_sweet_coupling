#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

if [ "x" != "x$DISPLAY" ]; then 
	echo
	echo "SPH VISUALIZATION"
	SCONS="scons --program=programs/vis_sphericalHarmonics --mode=debug"
	echo "$SCONS"
	$SCONS  || exit
fi


mule.benchmark.cleanup_all || exit 1
