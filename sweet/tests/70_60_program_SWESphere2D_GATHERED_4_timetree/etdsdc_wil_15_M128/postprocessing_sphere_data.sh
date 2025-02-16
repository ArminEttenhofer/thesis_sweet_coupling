#!/bin/bash

filename=job_bench_COMP_spspec_spdeal_fft_benchtime_thomp_release_RT_M2_LEGENDRE_N3_tsm_ETDSDC_lg_ADDT_lb_lc_n__tso4_tsob4_dt00120.00_M0128/output_prog_phi_pert_t00001296000.00000000.sweet
geometry=sphere
plot_type=map
extent_type=full
compare_solution=0
filename_ref=no_ref
roll=1
lim_ref=0
geopotential_elevation=1

./plot_sphere_data.py $filename $geometry $plot_type $extent_type $compare_solution $filename_ref $roll $lim_ref $geopotential_elevation
