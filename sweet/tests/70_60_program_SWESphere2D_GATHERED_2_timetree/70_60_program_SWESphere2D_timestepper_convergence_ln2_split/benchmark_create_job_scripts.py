#! /usr/bin/env python3

import sys
import copy
from mule.parHelper import *

from mule.JobMule import *
jg = JobGeneration()

from mule.rexi.REXICoefficients import *
from mule.rexi.trexi.TREXI import *
from mule.rexi.cirexi.CIREXI import *
from mule.rexi.brexi.BREXI import *



#
# Run simulation on cart2d or sphere2d
#
jg.compile.program = 'programs/PDE_SWESphere2D'

jg.unique_id_filter = ['runtime.simparams', 'parallelization', 'runtime.benchmark', 'runtime.rexi_params']


# Verbosity mode
jg.runtime.verbosity = 2

#
# Mode and Grid resolution
#
jg.runtime.space_res_spectral = 64
jg.runtime.space_res_grid = None

#jg.runtime.benchmark_name = "gaussian_bumps_pvd"
jg.runtime.benchmark_name = "gaussian_bump"

#
# Compute error
#
jg.runtime.compute_errors = 0

#
# Preallocate the REXI matrices
#
jg.runtime.rexi_sphere2d_preallocation = 1

#
# Threading across all REXI terms
#
jg.compile.rexi_thread_parallel_sum = "enable"



jg.runtime.f_sphere2d = 0

jg.runtime.viscosity = 0.0





ref_ts_size = 2
ref_ts_order = 4
ref_ts_method = f"ERK(ln,order={ref_ts_order})"

#
# A 2nd order accurate method already considerably reduces the errors
# Therefore, we use larger time step sizes to increase the errors
# to get errors larger than numerical precision
#
# We still want to have a very small time step size for the reference solution
# This is in particular important for REXI comparisons with ln2-type tests
#


ts_order = 2
o=f"order={ts_order}"
ts_methods = [
        f"ERK(ADDT(l,n),{o})",
        f"ERK(ADDT(lg,lc,n),{o})",

        #f"ln_erk_split_uv",
        f"ERK(ADDT(lg,lc,na_uv,nr_uv),{o})",
        #f"ln_erk_split_aa_uv",

        #f"ln_erk_split_vd",
        f"ERK(ADDT(lg,lc,na_vd,nr_vd),{o})",

        #f"ln_erk_split_aa_vd",
    ]


timestep_size_min = 64
timestep_sizes = [timestep_size_min*(2.0**i) for i in range(0, 6)]

jg.runtime.max_simulation_time = timestep_size_min*512

#####################################################
#####################################################
#####################################################

jg.runtime.output_timestep_size = jg.runtime.max_simulation_time



#
# Reference solution
#
jg.runtime.rexi_method = None
jg.runtime.timestepping_method = ref_ts_method
jg.runtime.timestepping_order = ref_ts_order
jg.runtime.timestepping_order2 = ref_ts_order
jg.runtime.timestep_size = ref_ts_size

jg.reference_job = True
jg.gen_jobscript_directory()
jg.reference_job = False

# Use this one as the reference solution!
jg.reference_job_unique_id = jg.job_unique_id



#
# Create job scripts
#

jg.runtime.timestepping_order = ts_order
jg.runtime.timestepping_order2 = ts_order


jg_base = copy.deepcopy(jg)

for tsm in ts_methods:
    for timestep_size in timestep_sizes:

        # Create new job description
        jg = copy.deepcopy(jg_base)

        jg.runtime.timestepping_method = tsm
        jg.runtime.timestep_size = timestep_size

        if jg.runtime.max_simulation_time % jg.runtime.timestep_size != 0:
            print("simtime: "+str(jg.runtime.max_simulation_time))
            print("timestep_size: "+str(jg.runtime.timestep_size))
            raise Exception("Invalid time step size (not remainder-less dividable)")

        if '_exp' in jg.runtime.timestepping_method or 'REXI(' in jg.runtime.timestepping_method:

            if 'rexi' in jg.runtime.timestepping_method.lower():
                # CI REXI via file input
                jg.runtime.rexi_method = 'file'

                cirexi = CIREXI()
                coeffs_phi0 = cirexi.setup("phi0", N=32, R=2).toFloat()
                coeffs_phi1 = cirexi.setup("phi1", N=32, R=2).toFloat()
                coeffs_phi2 = cirexi.setup("phi2", N=32, R=2).toFloat()
                coeffs_phi3 = cirexi.setup("phi3", N=32, R=2).toFloat()

                coeffs_phi0.normalize_steady_state()
                coeffs_phi1.normalize_steady_state()
                coeffs_phi2.normalize_steady_state()
                coeffs_phi3.normalize_steady_state()

                jg.runtime.rexi_files_coefficients = [coeffs_phi0, coeffs_phi1, coeffs_phi2, coeffs_phi3]

            else:
                # B REXI via file
                jg.runtime.rexi_method = 'file'

                brexi = BREXI()
                coeffs = brexi.setup(N=8, quadrature_method='gauss_legendre').toFloat()
                jg.runtime.rexi_files_coefficients = [coeffs]


        setupParallelization(jg)
        jg.gen_jobscript_directory()
