#! /usr/bin/env bash

precice=""

while getopts p: flag
do
    case "${flag}" in
        p) precice="-p ${OPTARG}";;
        *) ;;
    esac
done

rm -rf precice-run/
rm -rf output/
mkdir output

../sweet/build/programs/PDE_SWECart2D_COMP_plspec_pldeal_fft_thomp_release -N 200 --timestepping-method="ERK(ln, order=4)" --dt=1 --benchmark-name=column -t 1000 -o 5 -v 2 ${precice} --output-file-mode "csv" --output-file-name "output/output_%s_t%020.8f.csv" --instability-checks=true
