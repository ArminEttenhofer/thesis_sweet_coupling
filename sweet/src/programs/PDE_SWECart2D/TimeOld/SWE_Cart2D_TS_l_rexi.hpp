/*
 * SWE_Cart2D_TS_ln_erk.hpp
 *
 *  Created on: 29 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWECART2D_TIME_SWE_CART2D_TS_L_REXI_HPP
#define PROGRAMS_PDE_SWECART2D_TIME_SWE_CART2D_TS_L_REXI_HPP

#include <programs/PDE_SWECart2D/TimeOld/PDESWECart2DTS_BaseInterface.hpp>
#include <programs/PDE_SWECart2D/TimeOld/SWE_Cart2D_TS_l_direct.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Data/Cart2D/Operators.hpp>
#include <sweet/Data/Cart2DComplex/DataSpectral.hpp>
#include <sweet/Data/Cart2DComplex/Operators.hpp>
#include <limits>
#include <string>
#include <complex>
#include <sweet/ExpIntegration/Shack.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>



#if SWEET_BENCHMARK_TIMINGS
#	include <sweet/Tools/Stopwatch.hpp>
#endif


#if SWEET_MPI
#	include <mpi.h>
#endif

class SWE_Cart2D_TS_l_rexi	:
		public PDESWECart2DTS_BaseInterface
{
	std::vector<std::complex<double>> rexi_alphas;
	std::vector<std::complex<double>> rexi_betas;
	std::complex<double> rexi_gamma;

	//! simulation domain size
	double domain_size[2];

	std::size_t block_size;

#if SWEET_BENCHMARK_TIMINGS
	sweet::Tools::Stopwatch stopwatch_preprocessing;
	sweet::Tools::Stopwatch stopwatch_broadcast;
	sweet::Tools::Stopwatch stopwatch_reduce;
	sweet::Tools::Stopwatch stopwatch_solve_rexi_terms;
#endif

	class PerThreadVars
	{
	public:
		sweet::Data::Cart2DComplex::Operators ops;

		sweet::Data::Cart2DComplex::DataSpectral eta;

		sweet::Data::Cart2DComplex::DataSpectral eta0;
		sweet::Data::Cart2DComplex::DataSpectral u0;
		sweet::Data::Cart2DComplex::DataSpectral v0;

		sweet::Data::Cart2DComplex::DataSpectral h_sum;
		sweet::Data::Cart2DComplex::DataSpectral u_sum;
		sweet::Data::Cart2DComplex::DataSpectral v_sum;
	};

	//! per-thread allocated variables to avoid NUMA domain effects
	std::vector<PerThreadVars*> perThreadVars;

	//! number of threads to be used
	int num_local_rexi_par_threads;

	//! number of mpi ranks to be used
	int mpi_rank;

	//! MPI ranks
	int num_mpi_ranks;

	//! number of threads to be used
	int num_global_threads;

public:
	//! final time step
	bool final_timestep;

	//! use direct solution instead of REXI
	bool rexi_use_direct_solution;

	//! Direct solution for linear parts
	SWE_Cart2D_TS_l_direct ts_l_direct;

public:
	bool setup(
			sweet::Data::Cart2D::Operators *io_ops,
			const std::string &i_function_name
	);

	bool setup(
			sweet::Data::Cart2D::Operators *io_ops
	);

	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	);

	void runTimestep(
			const sweet::Data::Cart2D::DataSpectral &i_h_pert,	//!< prognostic variables
			const sweet::Data::Cart2D::DataSpectral &i_u,	//!< prognostic variables
			const sweet::Data::Cart2D::DataSpectral &i_v,	//!< prognostic variables

			sweet::Data::Cart2D::DataSpectral &o_h_pert,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &o_u,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &o_v,	//!< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);


	void run_timestep_real(
			const sweet::Data::Cart2D::DataSpectral &i_h_pert,	//!< prognostic variables
			const sweet::Data::Cart2D::DataSpectral &i_u,	//!< prognostic variables
			const sweet::Data::Cart2D::DataSpectral &i_v,	//!< prognostic variables

			sweet::Data::Cart2D::DataSpectral &o_h_pert,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &o_u,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &o_v,	//!< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);

	void runTimestep(
			sweet::Data::Cart2D::DataSpectral &io_h_pert,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);



	void cleanup();



public:
	static
	void MPI_quitWorkers(
			sweet::Data::Cart2D::Config *i_cart2DDataConfig
	);



	virtual ~SWE_Cart2D_TS_l_rexi();
};

#endif
