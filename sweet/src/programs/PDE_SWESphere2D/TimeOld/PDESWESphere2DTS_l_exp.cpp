/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Convert/DataSpectral_2_Sphere2DComplex_DataSpectral.hpp>
#include <sweet/Data/Sphere2DComplex/Convert/DataSpectral_2_Sphere2D_DataSpectral.hpp>
#include "PDESWESphere2DTS_l_exp.hpp"

#include <iostream>
#include <utility>
#include <sweet/Tools/StopwatchBox.hpp>
#include <sweet/ExpIntegration/ExpFunction.hpp>
#include <sweet/ExpIntegration/REXI/REXI.hpp>
#include "PDESWESphere2DTS_l_exp_direct_special.hpp"

#ifndef SWEET_THREADING_TIME_REXI
#	define SWEET_THREADING_TIME_REXI 1
#endif

#define SWEET_REXI_SPECTRAL_SPACE_REDUCTION	1

/*
 * Compute the REXI sum massively parallel *without* a parallelization with parfor in space
 */
#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
#	include <omp.h>
#endif


#ifndef SWEET_MPI
#	define SWEET_MPI 1
#endif


#if SWEET_MPI
#	include <mpi.h>
#endif


bool PDESWESphere2DTS_l_exp::shackRegistration(sweet::Shacks::Dictionary *io_shackDict)
{
	PDESWESphere2DTS_BaseInterface::shackRegistration(io_shackDict);

	timestepping_method_l_exp_direct_special.shackRegistration(shackDict);
	timestepping_method_lg_exp_lc_exp.shackRegistration(shackDict);

	return true;
}

bool PDESWESphere2DTS_l_exp::implementsTimesteppingMethod(const std::string &i_timestepping_method)
{
	timestepping_method = i_timestepping_method;
	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
	if (i_timestepping_method == "l_exp" || i_timestepping_method == "lg_exp")
		return true;

	return false;
}


std::string PDESWESphere2DTS_l_exp::getIDString()
{
	return "l_exp";
}


bool PDESWESphere2DTS_l_exp::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	bool no_coriolis = false;

	if (timestepping_method == "lg_exp")
		no_coriolis = true;

	return setup_variant_100(
		io_ops,
		shackExpIntegration,
		"phi0",
		shackExpIntegration->exp_method,
		shackTimestepControl->currentTimestepSize,
		shackPDESWESphere2D->sphere2d_use_fsphere2D,
		no_coriolis,
		shackPDESWETimeDisc->timestepping_order,
		shackExpIntegration->sphere2d_solver_preallocation
	);

}

bool PDESWESphere2DTS_l_exp::setup_variant_10(
		const sweet::Data::Sphere2D::Operators *io_ops,
		sweet::ExpIntegration::Shack *i_shackExpIntegration,
		const std::string &i_function_name,
		double i_timestepSize,
		bool i_use_f_sphere2D,
		bool i_no_coriolis
)
{
	return setup_variant_100(
		io_ops,
		shackExpIntegration,
		"phi0",
		shackExpIntegration->exp_method,
		i_timestepSize,
		i_use_f_sphere2D,
		i_no_coriolis,
		-1,
		shackExpIntegration->sphere2d_solver_preallocation
	);
}

bool PDESWESphere2DTS_l_exp::setup_variant_50(
		const sweet::Data::Sphere2D::Operators *io_ops,
		sweet::ExpIntegration::Shack *i_shackExpIntegration,
		const std::string &i_function_name,
		double i_timestepSize,
		bool i_use_f_sphere2D,
		bool i_no_coriolis,
		int i_timestepping_order
)
{
	return setup_variant_100(
		io_ops,
		i_shackExpIntegration,
		i_function_name,
		shackExpIntegration->exp_method,
		i_timestepSize,
		i_use_f_sphere2D,
		i_no_coriolis,
		i_timestepping_order,
		shackExpIntegration->sphere2d_solver_preallocation
	);
}


bool PDESWESphere2DTS_l_exp::setup_variant_100(
		const sweet::Data::Sphere2D::Operators *io_ops,
		sweet::ExpIntegration::Shack *i_shackExpIntegration,
		const std::string &i_functionName,
		const std::string &i_exp_method,
		double i_timestepSize,
		bool i_use_f_sphere2D,
		bool i_no_coriolis,
		int i_timestepping_order,
		bool i_use_rexi_sphere2d_solver_preallocation
)
{
	ops = io_ops;

	function_name = i_functionName;
	exp_method = i_exp_method;
	timestep_size = i_timestepSize;

	use_f_sphere2D = i_use_f_sphere2D;
	no_coriolis = i_no_coriolis;

	timestepping_order = i_timestepping_order;

	use_rexi_sphere2d_solver_preallocation = i_use_rexi_sphere2d_solver_preallocation;


	sweet::ExpIntegration::REXI::REXI<> rexi;

	/*
	 * Print some useful information
	 */
	if (	exp_method != "direct" &&
			exp_method != "ss_taylor"
	)
	{
		if (!rexi.is_rexi_method_supported(exp_method))
		{
			std::cerr << std::endl;
			std::cerr << "Available EXP methods (--exp-method=...):" << std::endl;
			std::cerr << "        'direct': Analytical solution" << std::endl;
			std::cerr << "        'ss_taylor': Strang-Split Taylor of exp(L) \approx exp(L_g) taylor(L_c)" << std::endl;

			std::stringstream ss;
			rexi.get_available_rexi_methods(ss);
			std::cerr << ss.str() << std::endl;

			if (exp_method == "help")
				SWEETErrorFatal("See above for available REXI methods");
			else
				SWEETErrorFatal("Unknown EXP method");
		}
	}


	#if SWEET_MPI
		if (mpi_rank == 0)
		{
			#if SWEET_REXI_SPECTRAL_SPACE_REDUCTION
				int rexi_communication_size = ops->sphere2DDataConfig->spectral_array_data_number_of_elements*2*sizeof(double);
			#else
				int rexi_communication_size = ops->sphere2DDataConfig->grid_number_elements*sizeof(double);
			#endif
			std::cout << "[MULE] rexi.communication_size: " << rexi_communication_size << std::endl;
		}
	#endif

	reset();

	#if SWEET_BENCHMARK_TIMINGS
		sweet::Tools::StopwatchBox::getInstance().rexi.start();
		sweet::Tools::StopwatchBox::getInstance().rexi_setup.start();
	#endif


	/*
	 * Setup REXI function evaluations
	 */
	_expFunction.setup(i_functionName);


	use_exp_method_direct_solution = false;
	use_exp_method_strang_split_taylor = false;
	use_exp_method_rexi = false;

	if (exp_method == "direct")
	{
		if (use_f_sphere2D)
			SWEETErrorFatal("f-sphere2D solution not implemented");

		use_exp_method_direct_solution = true;

		timestepping_method_l_exp_direct_special.setup_main(ops, timestepping_order, !no_coriolis, function_name);
	}
	else if (exp_method == "ss_taylor")
	{
		if (no_coriolis)
			SWEETErrorFatal("'ss_taylor' intended to include Coriolis effect!");

		use_exp_method_strang_split_taylor = true;

		timestepping_method_lg_exp_lc_exp.setup_main(ops, timestepping_order);
	}
	else
	{
		use_exp_method_rexi = true;


		sweet::ExpIntegration::REXI::REXICoefficients<double> rexiCoefficients;
		bool retval = rexi.load(
				i_shackExpIntegration,
				function_name,
				rexiCoefficients,
				shackExpIntegration->verbosity
		);

		if (!retval)
			SWEETErrorFatal(std::string("Phi function '")+function_name+std::string("' not provided or not supported"));

		_rexi_alphas = rexiCoefficients.alphas;
		for (std::size_t n = 0; n < _rexi_alphas.size(); n++)
			_rexi_alphas[n] = -_rexi_alphas[n];

		_rexi_betas = rexiCoefficients.betas;
		_rexi_gamma = rexiCoefficients.gamma;

		std::size_t N = _rexi_alphas.size();
		block_size = N/num_global_threads;
		if (block_size*num_global_threads != N)
			block_size++;

		perThreadVars.resize(num_local_rexi_par_threads);


		/**
		 * We split the setup from the utilization here.
		 *
		 * This is necessary, since it has to be assured that
		 * the FFTW plans are initialized before using them.
		 */
		if (num_local_rexi_par_threads == 0)
		{
			std::cerr << "FATAL ERROR B: omp_get_max_threads == 0" << std::endl;
			exit(-1);
		}

#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
		if (omp_in_parallel())
			SWEETErrorFatal("FATAL ERROR X: in parallel region");
#endif

		/*
		 * Use a kind of serialization of the input to avoid threading
		 * conflicts in the ComplexFFT generation
		 */
		for (int j = 0; j < num_local_rexi_par_threads; j++)
		{
			#if SWEET_THREADING_TIME_REXI
			#pragma omp parallel for schedule(static,1) default(none) shared(std::cout,std::cerr,j)
			#endif
			for (int local_thread_id = 0; local_thread_id < num_local_rexi_par_threads; local_thread_id++)
			{
				if (local_thread_id != j)
					continue;

				#if SWEET_DEBUG && SWEET_THREADING_TIME_REXI
					if (omp_get_thread_num() != local_thread_id)
					{
						// leave this dummy std::cout in it to avoid the intel compiler removing this part
						std::cout << "ERROR: thread " << omp_get_thread_num() << " number mismatch " << local_thread_id << std::endl;
						exit(-1);
					}
				#endif

				perThreadVars[local_thread_id] = new PerThreadVars;

				std::size_t start, end;
				p_get_workload_start_end(start, end, local_thread_id);
				int local_size = (int)end-(int)start;

				#if SWEET_DEBUG
					if (local_size < 0)
					{
						std::cerr << "local_size < 0" << std::endl;
						exit(-1);
					}
				#endif

				perThreadVars[local_thread_id]->alpha.resize(local_size);
				perThreadVars[local_thread_id]->beta.resize(local_size);

				perThreadVars[local_thread_id]->accum_phi.setup(ops->sphere2DDataConfig);
				perThreadVars[local_thread_id]->accum_vrt.setup(ops->sphere2DDataConfig);
				perThreadVars[local_thread_id]->accum_div.setup(ops->sphere2DDataConfig);

				for (std::size_t n = start; n < end; n++)
				{
					int thread_local_idx = n-start;

					perThreadVars[local_thread_id]->alpha[thread_local_idx] = _rexi_alphas[n];
					perThreadVars[local_thread_id]->beta[thread_local_idx] = _rexi_betas[n];
				}
			}
		}

		p_update_coefficients();

		if (num_local_rexi_par_threads == 0)
		{
			std::cerr << "FATAL ERROR C: omp_get_max_threads == 0" << std::endl;
			exit(-1);
		}

	}

	#if SWEET_BENCHMARK_TIMINGS
		sweet::Tools::StopwatchBox::getInstance().rexi_setup.stop();
		sweet::Tools::StopwatchBox::getInstance().rexi.stop();
	#endif


#if 0
	if (shackDict.misc.verbosity > 2)
	{
		if (shackDict.REXI.exp_method != "direct")
		{
			std::cout << "ALPHA:" << std::endl;
			for (std::size_t n = 0; n < rexi_alphas.size(); n++)
				std::cout << rexi_alphas[n] << std::endl;

			std::cout << "BETA:" << std::endl;
			for (std::size_t n = 0; n < _rexi_betas.size(); n++)
				std::cout << _rexi_betas[n] << std::endl;

			std::cout << "GAMMA:" << std::endl;
			std::cout << _rexi_gamma << std::endl;
		}
	}
#endif

	return true;
}




void PDESWESphere2DTS_l_exp::runTimestep(
	const sweet::Data::Sphere2D::DataSpectral &i_prog_phi0,
	const sweet::Data::Sphere2D::DataSpectral &i_prog_vrt0,
	const sweet::Data::Sphere2D::DataSpectral &i_prog_div0,

	sweet::Data::Sphere2D::DataSpectral &o_prog_phi0,
	sweet::Data::Sphere2D::DataSpectral &o_prog_vrt0,
	sweet::Data::Sphere2D::DataSpectral &o_prog_div0,

	double i_fixed_dt,
	double i_simulation_timestamp
)
{
	o_prog_phi0 = i_prog_phi0;
	o_prog_vrt0 = i_prog_vrt0;
	o_prog_div0 = i_prog_div0;

	runTimestep(o_prog_phi0, o_prog_vrt0, o_prog_div0, i_fixed_dt, i_simulation_timestamp);
}



PDESWESphere2DTS_l_exp::PDESWESphere2DTS_l_exp()	:
	use_rexi_sphere2d_solver_preallocation(false),
	use_exp_method_direct_solution(false),
	use_exp_method_strang_split_taylor(false),
	use_exp_method_rexi(false)
{
	#if SWEET_BENCHMARK_TIMINGS
		sweet::Tools::StopwatchBox::getInstance().rexi.start();
		sweet::Tools::StopwatchBox::getInstance().rexi_setup.start();
	#endif

	#if !SWEET_USE_LIBFFT
		SWEETErrorFatal("Spectral space required for solvers, use compile option --libfft=enable");
	#endif

	#if SWEET_THREADING_TIME_REXI

		num_local_rexi_par_threads = omp_get_max_threads();

		if (num_local_rexi_par_threads == 0)
		{
			std::cerr << "FATAL ERROR: omp_get_max_threads == 0" << std::endl;
			exit(-1);
		}
	#else
		num_local_rexi_par_threads = 1;
	#endif

	#if SWEET_MPI
		mpi_comm = MPI_COMM_WORLD;	// TODO: Make me more flexible in future versions
		MPI_Comm_rank(mpi_comm, &mpi_rank);
		MPI_Comm_size(mpi_comm, &num_mpi_ranks);

		num_global_threads = num_local_rexi_par_threads * num_mpi_ranks;

	#else

		num_global_threads = num_local_rexi_par_threads;

	#endif


	#if SWEET_BENCHMARK_TIMINGS
		sweet::Tools::StopwatchBox::getInstance().rexi_setup.stop();
		sweet::Tools::StopwatchBox::getInstance().rexi.stop();
	#endif
}



void PDESWESphere2DTS_l_exp::reset()
{
	#if SWEET_BENCHMARK_TIMINGS
		sweet::Tools::StopwatchBox::getInstance().rexi.start();
		sweet::Tools::StopwatchBox::getInstance().rexi_setup.start();
	#endif

	for (std::vector<PerThreadVars*>::iterator iter = perThreadVars.begin(); iter != perThreadVars.end(); ++iter)
	{
		PerThreadVars* p = *iter;
		delete p;
	}

	perThreadVars.resize(0);

	#if SWEET_BENCHMARK_TIMINGS
		sweet::Tools::StopwatchBox::getInstance().rexi_setup.stop();
		sweet::Tools::StopwatchBox::getInstance().rexi.stop();
	#endif
}



PDESWESphere2DTS_l_exp::~PDESWESphere2DTS_l_exp()
{
	#if SWEET_BENCHMARK_TIMINGS
		sweet::Tools::StopwatchBox::getInstance().rexi.start();
		sweet::Tools::StopwatchBox::getInstance().rexi_shutdown.start();
	#endif

	for (std::vector<PerThreadVars*>::iterator iter = perThreadVars.begin(); iter != perThreadVars.end(); ++iter)
	{
		PerThreadVars* p = *iter;
		delete p;
	}

	#if SWEET_BENCHMARK_TIMINGS
		#if SWEET_MPI

			int num_ranks;
			MPI_Comm_size(mpi_comm, &num_ranks);

			if (num_ranks > 1)
			{
				/*
				 * Send broadcast information from 2nd rank to 1st rank
				 *
				 * This is required in case of buffered broadcasts from the 1st rank
				 * which makes MPI_Bcast to return immediately.
				 */
				if (mpi_rank == 1)
				{
					double data = sweet::Tools::StopwatchBox::getInstance().rexi_timestepping_broadcast.time;
					MPI_Send(&data, sizeof(double), MPI_BYTE, 0, 0, mpi_comm);
				}

				if (mpi_rank == 0)
				{
					MPI_Status status;
					MPI_Recv(&sweet::Tools::StopwatchBox::getInstance().rexi_timestepping_broadcast.time, sizeof(double), MPI_BYTE, 1, 0, mpi_comm, &status);
				}

			}
		#endif

		sweet::Tools::StopwatchBox::getInstance().rexi_shutdown.stop();
		sweet::Tools::StopwatchBox::getInstance().rexi.stop();
	#endif
}



void PDESWESphere2DTS_l_exp::p_get_workload_start_end(
		std::size_t &o_start,
		std::size_t &o_end,
		int i_local_thread_id
)
{
	std::size_t max_N = _rexi_alphas.size();

	#if SWEET_THREADING_TIME_REXI || SWEET_MPI

		#if SWEET_MPI
			int global_thread_id = i_local_thread_id + num_local_rexi_par_threads*mpi_rank;
		#else
			int global_thread_id = i_local_thread_id;
		#endif

		SWEET_ASSERT(block_size >= 0);
		SWEET_ASSERT(global_thread_id >= 0);

		o_start = std::min(max_N, block_size*global_thread_id);
		o_end = std::min(max_N, o_start+block_size);

	#else

		o_start = 0;
		o_end = max_N;

	#endif
}




void PDESWESphere2DTS_l_exp::p_update_coefficients()
{
	#if SWEET_THREADING_TIME_REXI
	#pragma omp parallel for schedule(static,1) default(none) shared(std::cout)
	#endif
	for (int local_thread_id = 0; local_thread_id < num_local_rexi_par_threads; local_thread_id++)
	{
		std::size_t start, end;
		p_get_workload_start_end(start, end, local_thread_id);
		int local_size = (int)end-(int)start;

		if (use_rexi_sphere2d_solver_preallocation)
		{
			perThreadVars[local_thread_id]->sweREXISolvers.resize(local_size);

			for (std::size_t n = start; n < end; n++)
			{
				int thread_local_idx = n-start;

				perThreadVars[local_thread_id]->sweREXISolvers[thread_local_idx].setup_vectorinvariant_progphivortdiv(
						ops->sphere2DDataConfig,
						shackSphere2DDataOps,

						perThreadVars[local_thread_id]->alpha[thread_local_idx],
						perThreadVars[local_thread_id]->beta[thread_local_idx],
						shackSphere2DDataOps->sphere2d_radius,
						shackPDESWESphere2D->sphere2d_rotating_coriolis_omega,
						shackPDESWESphere2D->sphere2d_fsphere2d_f0,
						shackPDESWESphere2D->h0 * shackPDESWESphere2D->gravitation,
						timestep_size,
						use_f_sphere2D,
						no_coriolis
					);
			}
		}
	}
}




/**
 * Solve the REXI of \f$ U(t) = exp(L*t) \f$
 *
 * See
 * 		doc/rexi/understanding_rexi.pdf
 *
 * for further information
 */
void PDESWESphere2DTS_l_exp::runTimestep(
	sweet::Data::Sphere2D::DataSpectral &io_prog_phi,
	sweet::Data::Sphere2D::DataSpectral &io_prog_vrt,
	sweet::Data::Sphere2D::DataSpectral &io_prog_div,

	double i_fixed_dt,
	double i_simulation_timestamp
)
{
	if (use_exp_method_direct_solution)
	{
		#if SWEET_BENCHMARK_TIMINGS
			sweet::Tools::StopwatchBox::getInstance().rexi.start();
			sweet::Tools::StopwatchBox::getInstance().rexi_timestepping.start();
		#endif

		timestepping_method_l_exp_direct_special.runTimestep(io_prog_phi, io_prog_vrt, io_prog_div, i_fixed_dt, i_simulation_timestamp);
#if 0
		// no Coriolis force active

		/*
		 * Using exponential integrators, we must compute an
		 */

		double ir = 1.0/shackSphere2DDataOps->sphere2d_radius;
		// avg. geopotential

		double G = -simCoeffs.h0*simCoeffs.gravitation;

		/*
		 * See doc/rexi/rexi_for_swe_on_nonrotating_sphere2D.pdf
		 */

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			std::size_t idx = sphere2DDataConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				double D = -(double)n*((double)n+1.0)*ir*ir;
				D = -D;

				if (D == 0)
				{
					idx++;
					continue;
				}

				std::complex<double> &phi0 = io_prog_phi.spectral_space_data[idx];
				std::complex<double> &div0 = io_prog_div.spectral_space_data[idx];

				// TODO: precompute this

				// result will be imaginary only!
				std::complex<double> sqrt_DG = std::sqrt(std::complex<double>(D*G));

				// Multiply with Q^{-1}
				std::complex<double> l0 = -sqrt_DG/(2*G) * phi0 + 0.5*div0;
				std::complex<double> l1 = +sqrt_DG/(2*G) * phi0 + 0.5*div0;

				l0 = _expFunction.eval(timestep_size*(-sqrt_DG))*l0;
				l1 = _expFunction.eval(timestep_size*sqrt_DG)*l1;

				phi0 = -G/sqrt_DG * l0 + G/sqrt_DG* l1;
				div0 = l0 + l1;

				idx++;
			}
		}
#endif


		#if SWEET_BENCHMARK_TIMINGS
			sweet::Tools::StopwatchBox::getInstance().rexi_timestepping.stop();
			sweet::Tools::StopwatchBox::getInstance().rexi.stop();
		#endif

		return;
	}

	if (use_exp_method_strang_split_taylor)
	{
		#if SWEET_BENCHMARK_TIMINGS
			sweet::Tools::StopwatchBox::getInstance().rexi.start();
			sweet::Tools::StopwatchBox::getInstance().rexi_timestepping.start();
		#endif

		timestepping_method_lg_exp_lc_exp.runTimestep(io_prog_phi, io_prog_vrt, io_prog_div, i_fixed_dt, i_simulation_timestamp);

		#if SWEET_BENCHMARK_TIMINGS
			sweet::Tools::StopwatchBox::getInstance().rexi_timestepping.stop();
			sweet::Tools::StopwatchBox::getInstance().rexi.stop();
		#endif

		return;
	}

	if (use_exp_method_strang_split_taylor)
	{
		#if SWEET_BENCHMARK_TIMINGS
			sweet::Tools::StopwatchBox::getInstance().rexi.start();
			sweet::Tools::StopwatchBox::getInstance().rexi_timestepping.start();
		#endif

		timestepping_method_l_exp_direct_special.runTimestep(io_prog_phi, io_prog_vrt, io_prog_div, i_fixed_dt, i_simulation_timestamp);
#if 0
		// no Coriolis force active

		/*
		 * Using exponential integrators, we must compute an
		 */

		double ir = 1.0/shackSphere2DDataOps->sphere2d_radius;
		// avg. geopotential

		double G = -simCoeffs.h0*simCoeffs.gravitation;

		/*
		 * See doc/rexi/rexi_for_swe_on_nonrotating_sphere2D.pdf
		 */

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			std::size_t idx = sphere2DDataConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				double D = -(double)n*((double)n+1.0)*ir*ir;
				D = -D;

				if (D == 0)
				{
					idx++;
					continue;
				}

				std::complex<double> &phi0 = io_prog_phi.spectral_space_data[idx];
				std::complex<double> &div0 = io_prog_div.spectral_space_data[idx];

				// TODO: precompute this

				// result will be imaginary only!
				std::complex<double> sqrt_DG = std::sqrt(std::complex<double>(D*G));

				// Multiply with Q^{-1}
				std::complex<double> l0 = -sqrt_DG/(2*G) * phi0 + 0.5*div0;
				std::complex<double> l1 = +sqrt_DG/(2*G) * phi0 + 0.5*div0;

				l0 = _expFunction.eval(timestep_size*(-sqrt_DG))*l0;
				l1 = _expFunction.eval(timestep_size*sqrt_DG)*l1;

				phi0 = -G/sqrt_DG * l0 + G/sqrt_DG* l1;
				div0 = l0 + l1;

				idx++;
			}
		}
#endif


		#if SWEET_BENCHMARK_TIMINGS
			sweet::Tools::StopwatchBox::getInstance().rexi_timestepping.stop();
			sweet::Tools::StopwatchBox::getInstance().rexi.stop();
		#endif

		return;
	}


#if SWEET_BENCHMARK_TIMINGS
	sweet::Tools::StopwatchBox::getInstance().rexi.start();
	sweet::Tools::StopwatchBox::getInstance().rexi_timestepping.start();
#endif

	bool use_rexi_gamma = false;

#if SWEET_MPI
	if (mpi_rank == 0)
#endif
	{
		if (_rexi_gamma.real() != 0)
		{
			use_rexi_gamma = true;
		}
	}

	/*
	 * PREPROCESSING
	 */
	#if SWEET_BENCHMARK_TIMINGS
		sweet::Tools::StopwatchBox::getInstance().rexi_timestepping_miscprocessing.start();
	#endif

		if (i_fixed_dt <= 0)
			SWEETErrorFatal("Only constant time step size allowed");

		double update_dt_delta = std::abs(timestep_size - i_fixed_dt)/std::max(timestep_size, i_fixed_dt);
		if (update_dt_delta > 1e-9)
		{
			std::cout << "Warning: Updating time step size from " << i_fixed_dt << " to " << timestep_size << ", threshold " << update_dt_delta << " exceeded" << std::endl;
			std::cout << "Warning: This warning is only triggered because certain things need to be recomputed" << std::endl;

			std::cout << "timestep_size: " << timestep_size << std::endl;
			std::cout << "i_fixed_dt: " << i_fixed_dt << std::endl;
			std::cout << "a: " << std::abs(timestep_size - i_fixed_dt) << std::endl;
			std::cout << "b: " << std::max(timestep_size, i_fixed_dt) << std::endl;
			std::cout << "c: " << std::abs(timestep_size - i_fixed_dt)/std::max(timestep_size, i_fixed_dt) << std::endl;
			std::cout << "update_dt_delta: " << update_dt_delta << std::endl;

			timestep_size = i_fixed_dt;

			//p_update_coefficients(true);
			p_update_coefficients();
		}

		#if SWEET_REXI_TIMINGS_ADDITIONAL_BARRIERS && SWEET_MPI
			MPI_Barrier(mpi_comm);
		#endif

	#if SWEET_BENCHMARK_TIMINGS
		sweet::Tools::StopwatchBox::getInstance().rexi_timestepping_miscprocessing.stop();
	#endif

	/*
	 * BROADCAST
	 */

	#if SWEET_BENCHMARK_TIMINGS
		sweet::Tools::StopwatchBox::getInstance().rexi_timestepping_broadcast.start();
	#endif

		#if SWEET_MPI
			/*
			 * We should measure this for the 2nd rank! And we do so (see later on)
			 */

			std::size_t spectral_data_num_doubles = io_prog_phi.sphere2DDataConfig->spectral_array_data_number_of_elements*2;

			MPI_Bcast(io_prog_phi.spectral_space_data, spectral_data_num_doubles, MPI_DOUBLE, 0, mpi_comm);
			MPI_Bcast(io_prog_vrt.spectral_space_data, spectral_data_num_doubles, MPI_DOUBLE, 0, mpi_comm);
			MPI_Bcast(io_prog_div.spectral_space_data, spectral_data_num_doubles, MPI_DOUBLE, 0, mpi_comm);

		#endif

	#if SWEET_BENCHMARK_TIMINGS
		sweet::Tools::StopwatchBox::getInstance().rexi_timestepping_broadcast.stop();
	#endif


	/*
	 * Special handler for REXI without using mode extension and no threading to speedup things
	 */
	#if !SWEET_THREADING_TIME_REXI

		#if SWEET_BENCHMARK_TIMINGS
			sweet::Tools::StopwatchBox::getInstance().rexi_timestepping_solver.start();
		#endif

			std::size_t start, end;
			p_get_workload_start_end(start, end, 0);

			perThreadVars[0]->accum_phi.spectral_setZero();
			perThreadVars[0]->accum_vrt.spectral_setZero();
			perThreadVars[0]->accum_div.spectral_setZero();

			{
				/*
				 * -> No threading
				 */

				sweet::Data::Sphere2D::DataSpectral tmp_prog_phi(ops->sphere2DDataConfig);
				sweet::Data::Sphere2D::DataSpectral tmp_prog_vort(ops->sphere2DDataConfig);
				sweet::Data::Sphere2D::DataSpectral tmp_prog_div(ops->sphere2DDataConfig);

				for (std::size_t workload_idx = start; workload_idx < end; workload_idx++)
				{
					int local_idx = workload_idx-start;

					if (use_rexi_sphere2d_solver_preallocation)
					{
						perThreadVars[0]->sweREXISolvers[local_idx].solve_vectorinvariant_progphivortdiv(
								io_prog_phi, io_prog_vrt, io_prog_div,
								tmp_prog_phi, tmp_prog_vort, tmp_prog_div
							);
					}
					else
					{
						SWERexiTerm_SPH rexiSPH;

						std::complex<double> &alpha = perThreadVars[0]->alpha[local_idx];
						std::complex<double> &beta = perThreadVars[0]->beta[local_idx];

						rexiSPH.setup_vectorinvariant_progphivortdiv(
								ops->sphere2DDataConfig,
								shackSphere2DDataOps,

								alpha,
								beta,

								shackSphere2DDataOps->sphere2d_radius,
								shackPDESWESphere2D->sphere2d_rotating_coriolis_omega,
								shackPDESWESphere2D->sphere2d_fsphere2d_f0,
								shackPDESWESphere2D->h0 * shackPDESWESphere2D->gravitation,
								i_fixed_dt,

								use_f_sphere2D,
								no_coriolis
						);

						rexiSPH.solve_vectorinvariant_progphivortdiv(
								io_prog_phi, io_prog_vrt, io_prog_div,
								tmp_prog_phi, tmp_prog_vort, tmp_prog_div
							);
					}

					perThreadVars[0]->accum_phi += tmp_prog_phi;
					perThreadVars[0]->accum_vrt += tmp_prog_vort;
					perThreadVars[0]->accum_div += tmp_prog_div;
				}


				if (use_rexi_gamma)
				{
					io_prog_phi = perThreadVars[0]->accum_phi + io_prog_phi*_rexi_gamma.real();
					io_prog_vrt = perThreadVars[0]->accum_vrt + io_prog_vrt*_rexi_gamma.real();
					io_prog_div = perThreadVars[0]->accum_div + io_prog_div*_rexi_gamma.real();
				}
				else
				{
					io_prog_phi = perThreadVars[0]->accum_phi;
					io_prog_vrt = perThreadVars[0]->accum_vrt;
					io_prog_div = perThreadVars[0]->accum_div;
				}

			}

		#if SWEET_BENCHMARK_TIMINGS
			sweet::Tools::StopwatchBox::getInstance().rexi_timestepping_solver.stop();
		#endif


	#else	/* SWEET_THREADING_TIME_REXI */


		{
			#if SWEET_BENCHMARK_TIMINGS
				sweet::Tools::StopwatchBox::getInstance().rexi_timestepping_solver.start();
			#endif

				#pragma omp parallel for schedule(static,1) default(none) shared(i_fixed_dt, io_prog_phi, io_prog_vrt, io_prog_div, std::cout, std::cerr)
				for (int local_thread_id = 0; local_thread_id < num_local_rexi_par_threads; local_thread_id++)
				{
					std::size_t start, end;
					p_get_workload_start_end(start, end, local_thread_id);

					/*
					* Make a copy to ensure that there are no race conditions by converting to physical space
					*/
					sweet::Data::Sphere2D::DataSpectral thread_io_prog_phi0 = io_prog_phi;
					sweet::Data::Sphere2D::DataSpectral thread_io_prog_vrt0 = io_prog_vrt;
					sweet::Data::Sphere2D::DataSpectral thread_io_prog_div0 = io_prog_div;

					sweet::Data::Sphere2D::DataSpectral tmp_prog_phi(ops->sphere2DDataConfig);
					sweet::Data::Sphere2D::DataSpectral tmp_prog_vort(ops->sphere2DDataConfig);
					sweet::Data::Sphere2D::DataSpectral tmp_prog_div(ops->sphere2DDataConfig);

					perThreadVars[local_thread_id]->accum_phi.spectral_setZero();
					perThreadVars[local_thread_id]->accum_vrt.spectral_setZero();
					perThreadVars[local_thread_id]->accum_div.spectral_setZero();


					for (std::size_t workload_idx = start; workload_idx < end; workload_idx++)
					{
						int local_idx = workload_idx-start;

						if (use_rexi_sphere2d_solver_preallocation)
						{
							perThreadVars[local_thread_id]->sweREXISolvers[local_idx].solve_vectorinvariant_progphivortdiv(
									thread_io_prog_phi0, thread_io_prog_vrt0, thread_io_prog_div0,
									tmp_prog_phi, tmp_prog_vort, tmp_prog_div
								);
						}
						else
						{
							SWERexiTerm_SPH rexiSPH;

							const std::complex<double> &alpha = perThreadVars[local_thread_id]->alpha[local_idx];
							const std::complex<double> &beta = perThreadVars[local_thread_id]->beta[local_idx];
							rexiSPH.setup_vectorinvariant_progphivortdiv(
									ops->sphere2DDataConfig,	//!< sphere2D data for input data
									shackSphere2DDataOps,

									alpha,
									beta,

									shackSphere2DDataOps->sphere2d_radius,
									shackPDESWESphere2D->sphere2d_rotating_coriolis_omega,
									shackPDESWESphere2D->sphere2d_fsphere2d_f0,
									shackPDESWESphere2D->h0*shackPDESWESphere2D->gravitation,
									i_fixed_dt,

									use_f_sphere2D,
									no_coriolis
							);

							rexiSPH.solve_vectorinvariant_progphivortdiv(
									io_prog_phi, io_prog_vrt, io_prog_div,
									tmp_prog_phi, tmp_prog_vort, tmp_prog_div
								);
						}

						perThreadVars[local_thread_id]->accum_phi += tmp_prog_phi;
						perThreadVars[local_thread_id]->accum_vrt += tmp_prog_vort;
						perThreadVars[local_thread_id]->accum_div += tmp_prog_div;
					}
				}

			#if SWEET_BENCHMARK_TIMINGS
				sweet::Tools::StopwatchBox::getInstance().rexi_timestepping_solver.stop();
				sweet::Tools::StopwatchBox::getInstance().rexi_timestepping_reduce.start();
			#endif

				if (use_rexi_gamma)
				{
					io_prog_phi *= _rexi_gamma.real();
					io_prog_vrt *= _rexi_gamma.real();
					io_prog_div *= _rexi_gamma.real();
				}
				else
				{
					io_prog_phi.spectral_setZero();
					io_prog_vrt.spectral_setZero();
					io_prog_div.spectral_setZero();
				}

				for (int thread_id = 0; thread_id < num_local_rexi_par_threads; thread_id++)
				{
					SWEET_ASSERT(io_prog_phi.sphere2DDataConfig->spectral_array_data_number_of_elements == perThreadVars[0]->accum_phi.sphere2DDataConfig->spectral_array_data_number_of_elements);

#if 0
					perThreadVars[thread_id]->accum_phi.request_data_physical();
					#pragma omp parallel for schedule(static) default(none) shared(io_prog_phi, thread_id)
					for (int i = 0; i < io_prog_phi.sphere2DDataConfig->grid_number_elements; i++)
						io_prog_phi.grid_space_data[i] += perThreadVars[thread_id]->accum_phi.grid_space_data[i];

					perThreadVars[thread_id]->accum_vrt.request_data_physical();
					#pragma omp parallel for schedule(static) default(none) shared(io_prog_vrt, thread_id)
					for (int i = 0; i < io_prog_vrt.sphere2DDataConfig->grid_number_elements; i++)
						io_prog_vrt.grid_space_data[i] += perThreadVars[thread_id]->accum_vrt.grid_space_data[i];


					perThreadVars[thread_id]->accum_div.request_data_physical();
					#pragma omp parallel for schedule(static) default(none) shared(io_prog_div, thread_id)
					for (int i = 0; i < io_prog_div0.sphere2DDataConfig->grid_number_elements; i++)
						io_prog_div.grid_space_data[i] += perThreadVars[thread_id]->accum_div.grid_space_data[i];
#else
					#pragma omp parallel for schedule(static) default(none) shared(io_prog_phi, thread_id)
					for (int i = 0; i < io_prog_phi.sphere2DDataConfig->spectral_array_data_number_of_elements; i++)
						io_prog_phi.spectral_space_data[i] += perThreadVars[thread_id]->accum_phi.spectral_space_data[i];

					#pragma omp parallel for schedule(static) default(none) shared(io_prog_vrt, thread_id)
					for (int i = 0; i < io_prog_vrt.sphere2DDataConfig->spectral_array_data_number_of_elements; i++)
						io_prog_vrt.spectral_space_data[i] += perThreadVars[thread_id]->accum_vrt.spectral_space_data[i];

					#pragma omp parallel for schedule(static) default(none) shared(io_prog_div, thread_id)
					for (int i = 0; i < io_prog_div.sphere2DDataConfig->spectral_array_data_number_of_elements; i++)
						io_prog_div.spectral_space_data[i] += perThreadVars[thread_id]->accum_div.spectral_space_data[i];
#endif
				}

			#if SWEET_BENCHMARK_TIMINGS
				sweet::Tools::StopwatchBox::getInstance().rexi_timestepping_reduce.stop();
			#endif
		}

	#endif	// END SWEET_THREADING_TIME_REXI

	#if SWEET_REXI_TIMINGS_ADDITIONAL_BARRIERS && SWEET_MPI
		#if SWEET_BENCHMARK_TIMINGS
			sweet::Tools::StopwatchBox::getInstance().rexi_timestepping_miscprocessing.start();
		#endif

			MPI_Barrier(mpi_comm);

		#if SWEET_BENCHMARK_TIMINGS
			sweet::Tools::StopwatchBox::getInstance().rexi_timestepping_miscprocessing.stop();
		#endif
	#endif


	#if SWEET_MPI

		#if SWEET_BENCHMARK_TIMINGS
			sweet::Tools::StopwatchBox::getInstance().rexi_timestepping_reduce.start();
		#endif

		#if SWEET_REXI_SPECTRAL_SPACE_REDUCTION

			/*
			 * Reduction in spectral space
			 *
			 * There was once a reason for this to be done in physical space.
			 * However, eventually, it seems to work also in spectral space.
			 */

			sweet::Data::Sphere2D::DataSpectral tmp(ops->sphere2DDataConfig);

// already setup
//			std::size_t spectral_data_num_doubles = io_prog_phi0.sphere2DDataConfig->spectral_array_data_number_of_elements*2;

			#if SWEET_REXI_ALLREDUCE

				MPI_Allreduce(prog_phi0_phys.spectral_space_data, tmp.spectral_space_data, spectral_data_num_doubles, MPI_DOUBLE, MPI_SUM, mpi_comm);
				std::swap(prog_phi0_phys.spectral_space_data, tmp.spectral_space_data);

				MPI_Allreduce(prog_vrt0_phys.spectral_space_data, tmp.spectral_space_data, spectral_data_num_doubles, MPI_DOUBLE, MPI_SUM, mpi_comm);
				std::swap(prog_vrt0.spectral_space_data, tmp.spectral_space_data);

				MPI_Allreduce(prog_div0_phys.spectral_space_data, tmp.spectral_space_data, spectral_data_num_doubles, MPI_DOUBLE, MPI_SUM, mpi_comm);
				std::swap(prog_div0.spectral_space_data, tmp.spectral_space_data);

			#else

				MPI_Reduce(io_prog_phi.spectral_space_data, tmp.spectral_space_data, spectral_data_num_doubles, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);
				if (mpi_rank == 0)
					std::swap(io_prog_phi.spectral_space_data, tmp.spectral_space_data);

				MPI_Reduce(io_prog_vrt.spectral_space_data, tmp.spectral_space_data, spectral_data_num_doubles, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);
				if (mpi_rank == 0)
					std::swap(io_prog_vrt.spectral_space_data, tmp.spectral_space_data);

				MPI_Reduce(io_prog_div.spectral_space_data, tmp.spectral_space_data, spectral_data_num_doubles, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);
				if (mpi_rank == 0)
					std::swap(io_prog_div.spectral_space_data, tmp.spectral_space_data);

			#endif

		#else

			/*
			 * Reduction in physical space
			 */

			sweet::Data::Sphere2D::DataGrid prog_phi0_phys = io_prog_phi.toGrid();
			sweet::Data::Sphere2D::DataGrid prog_vrt0_phys = io_prog_vrt.toGrid();
			sweet::Data::Sphere2D::DataGrid prog_div0_phys = io_prog_div.toGrid();

			/*
			 * Grid data reduction
			 *
			 * WE HAVE to do the reduction in physical space!
			 *
			 * Comment from Martin to Martin: I forgot why this was necessary :-(
			 */

			sweet::Data::Sphere2D::DataGrid tmp(sphere2DDataConfig);

			std::size_t grid_data_num_doubles = prog_phi0_phys.sphere2DDataConfig->grid_number_elements;

			#if SWEET_REXI_ALLREDUCE

				MPI_Allreduce(prog_phi0_phys.grid_space_data, tmp.grid_space_data, grid_data_num_doubles, MPI_DOUBLE, MPI_SUM, mpi_comm);
				std::swap(prog_phi0_phys.grid_space_data, tmp.grid_space_data);

				MPI_Allreduce(prog_vrt0_phys.grid_space_data, tmp.grid_space_data, grid_data_num_doubles, MPI_DOUBLE, MPI_SUM, mpi_comm);
				std::swap(prog_vrt0.grid_space_data, tmp.grid_space_data);

				MPI_Allreduce(prog_div0_phys.grid_space_data, tmp.grid_space_data, grid_data_num_doubles, MPI_DOUBLE, MPI_SUM, mpi_comm);
				std::swap(prog_div0.grid_space_data, tmp.grid_space_data);

			#else

				MPI_Reduce(prog_phi0_phys.grid_space_data, tmp.grid_space_data, grid_data_num_doubles, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);
				if (mpi_rank == 0)
					std::swap(prog_phi0_phys.grid_space_data, tmp.grid_space_data);

				MPI_Reduce(prog_vrt0_phys.grid_space_data, tmp.grid_space_data, grid_data_num_doubles, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);
				if (mpi_rank == 0)
					std::swap(prog_vrt0_phys.grid_space_data, tmp.grid_space_data);

				MPI_Reduce(prog_div0_phys.grid_space_data, tmp.grid_space_data, grid_data_num_doubles, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);
				if (mpi_rank == 0)
					std::swap(prog_div0_phys.grid_space_data, tmp.grid_space_data);

			#endif

			io_prog_phi.loadSphere2DDataGrid(prog_phi0_phys);
			io_prog_vrt.loadSphere2DDataGrid(prog_vrt0_phys);
			io_prog_div.loadSphere2DDataGrid(prog_div0_phys);

		#endif

		#if SWEET_BENCHMARK_TIMINGS
			sweet::Tools::StopwatchBox::getInstance().rexi_timestepping_reduce.stop();
		#endif
	#endif


	#if SWEET_BENCHMARK_TIMINGS
		sweet::Tools::StopwatchBox::getInstance().rexi_timestepping.stop();
		sweet::Tools::StopwatchBox::getInstance().rexi.stop();
	#endif
}

