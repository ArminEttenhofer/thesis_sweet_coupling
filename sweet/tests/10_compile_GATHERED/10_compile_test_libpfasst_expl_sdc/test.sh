#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

echo_info_hline
echo_info "PFASST-EXPL-SDC"
echo_info_hline


if [ -e "$MULE_LOCAL_ROOT/../local_software/local/lib/libpfasst.a" ]; then 

	SCONS="scons --program=programs/libpfasst/PDE_SWESphere2D_expl_sdc --libpfasst=enable --sweet-mpi=enable --libsph=enable --cart2d-spectral-space=disable --sphere2d-spectral-space=enable --threading=off --libfft=enable --sphere2d-spectral-dealiasing=enable"
	echo "$SCONS"
	$SCONS || exit

fi

mule.benchmark.cleanup_all || exit 1
