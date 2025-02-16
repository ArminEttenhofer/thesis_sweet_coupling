/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <programs/PDE_SWECart2D/TimeOld/SWE_Cart2D_TS_l_rexi.hpp>
///#include <sweet/Data/Cart2D/Convert/Complex_DataSpectral_2_DataSpectral.hpp>
///#include <sweet/Data/Cart2D/Convert/DataSpectral_2_Complex_DataSpectral.hpp>
#include <sweet/Data/Cart2D/Convert/DataSpectral_2_Cart2DComplex_DataSpectral.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Data/Cart2DComplex/Convert/DataSpectral_2_Cart2D_DataSpectral.hpp>
#include <sweet/Data/Cart2DComplex/DataSpectral.hpp>
#include <sweet/Data/Cart2DComplex/Operators.hpp>
#include <sweet/ExpIntegration/REXI/REXI.hpp>
#include <sweet/Tools/StopwatchBox.hpp>

#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
	#include <omp.h>
#endif


bool SWE_Cart2D_TS_l_rexi::setup(
	sweet::Data::Cart2D::Operators *io_ops
)
{
	return setup(io_ops, "phi0");
}


bool SWE_Cart2D_TS_l_rexi::setup(
	sweet::Data::Cart2D::Operators *io_ops,
	const std::string &i_function_name
)
{
	PDESWECart2DTS_BaseInterface::setup(io_ops);

	if (shackCart2DDataOps->space_grid_use_c_staggering)
		SWEETErrorFatal("Staggering not supported for l_rexi");

#if !SWEET_USE_LIBFFT
	std::cerr << "Spectral space required for solvers, use compile option --libfft=enable" << std::endl;
	exit(-1);
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
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_mpi_ranks);
#else
	mpi_rank = 0;
	num_mpi_ranks = 1;
#endif

	num_global_threads = num_local_rexi_par_threads * num_mpi_ranks;


	PDESWECart2DTS_BaseInterface::setup(io_ops);

	SWEET_ASSERT(shackCart2DDataOps != nullptr);
	SWEET_ASSERT(shackExpIntegration != nullptr);

	domain_size[0] = shackCart2DDataOps->cart2d_domain_size[0];
	domain_size[1] = shackCart2DDataOps->cart2d_domain_size[1];

	rexi_use_direct_solution = (shackExpIntegration->exp_method == "direct");

	if (rexi_use_direct_solution)
	{
		ts_l_direct.setup(io_ops, i_function_name);
		return true;
	}
	sweet::ExpIntegration::REXI::REXICoefficients<double> rexiCoefficients;


	sweet::ExpIntegration::REXI::REXI<> rexi;
	bool retval = rexi.load(
			shackExpIntegration,
			i_function_name,
			rexiCoefficients,
			shackIOData->verbosity
	);

	rexi_alphas = rexiCoefficients.alphas;
	rexi_betas = rexiCoefficients.betas;
	rexi_gamma = rexiCoefficients.gamma;


	if (!retval)
		SWEETErrorFatal(std::string("Phi function '")+i_function_name+std::string("' not available"));


	std::cout << "Number of total REXI coefficients N = " << rexi_alphas.size() << std::endl;

	std::size_t N = rexi_alphas.size();
	block_size = N/num_global_threads;
	if (block_size*num_global_threads != N)
		block_size++;

	cleanup();

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
	{
		std::cerr << "FATAL ERROR X: in parallel region" << std::endl;
		exit(-1);
	}
#endif

	sweet::Data::Cart2D::Config *cart2DDataConfig_local = ops->cart2DDataConfig;

	// use a kind of serialization of the input to avoid threading conflicts in the ComplexFFT generation
	for (int j = 0; j < num_local_rexi_par_threads; j++)
	{
#if SWEET_THREADING_TIME_REXI
#	pragma omp parallel for schedule(static,1) default(none) shared(cart2DDataConfig_local, std::cout,j)
#endif
		for (int i = 0; i < num_local_rexi_par_threads; i++)
		{
			if (i != j)
				continue;

#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
			if (omp_get_thread_num() != i)
			{
				// leave this dummy std::cout in it to avoid the intel compiler removing this part
				std::cout << "ERROR: thread " << omp_get_thread_num() << " number mismatch " << i << std::endl;
				exit(-1);
			}
#endif

			perThreadVars[i] = new PerThreadVars;

			perThreadVars[i]->ops.setup(cart2DDataConfig_local, domain_size);

			perThreadVars[i]->eta.setup(cart2DDataConfig_local);
			perThreadVars[i]->eta0.setup(cart2DDataConfig_local);
			perThreadVars[i]->u0.setup(cart2DDataConfig_local);
			perThreadVars[i]->v0.setup(cart2DDataConfig_local);
			perThreadVars[i]->h_sum.setup(cart2DDataConfig_local);
			perThreadVars[i]->u_sum.setup(cart2DDataConfig_local);
			perThreadVars[i]->v_sum.setup(cart2DDataConfig_local);
		}
	}

	if (num_local_rexi_par_threads == 0)
	{
		std::cerr << "FATAL ERROR C: omp_get_max_threads == 0" << std::endl;
		exit(-1);
	}

	for (int i = 0; i < num_local_rexi_par_threads; i++)
	{
		if (perThreadVars[i]->ops.diff_c_x.spectral_space_data == nullptr)
		{
			std::cerr << "ARRAY NOT INITIALIZED!!!!" << std::endl;
			exit(-1);
		}
	}

#if SWEET_THREADING_TIME_REXI
#	pragma omp parallel for schedule(static,1) default(none)  shared(cart2DDataConfig_local, std::cout)
#endif
	for (int i = 0; i < num_local_rexi_par_threads; i++)
	{
#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
		if (omp_get_thread_num() != i)
		{
			// leave this dummy std::cout in it to avoid the intel compiler removing this part
			std::cout << "ERROR: thread " << omp_get_thread_num() << " number mismatch " << i << std::endl;
			exit(-1);
		}
#else

#endif

		if (perThreadVars[i]->eta.spectral_space_data == nullptr)
		{
			std::cout << "ERROR, data == nullptr" << std::endl;
			exit(-1);
		}

		perThreadVars[i]->ops.setup(cart2DDataConfig_local, domain_size);

		// initialize all values to account for first touch policy reason
		perThreadVars[i]->eta.spectral_setZero();
		perThreadVars[i]->eta0.spectral_setZero();

		perThreadVars[i]->u0.spectral_setZero();
		perThreadVars[i]->v0.spectral_setZero();

		perThreadVars[i]->h_sum.spectral_setZero();
		perThreadVars[i]->u_sum.spectral_setZero();
		perThreadVars[i]->v_sum.spectral_setZero();

	}


#if SWEET_BENCHMARK_TIMINGS
	stopwatch_preprocessing.reset();
	stopwatch_broadcast.reset();
	stopwatch_reduce.reset();
	stopwatch_solve_rexi_terms.reset();
#endif

	return true;
}



void SWE_Cart2D_TS_l_rexi::runTimestep(
		const sweet::Data::Cart2D::DataSpectral &i_h_pert,	//!< prognostic variables
		const sweet::Data::Cart2D::DataSpectral &i_u,	//!< prognostic variables
		const sweet::Data::Cart2D::DataSpectral &i_v,	//!< prognostic variables

		sweet::Data::Cart2D::DataSpectral &o_h_pert,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &o_u,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &o_v,	//!< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	//! WARNING: i_h_pert might be identical to o_h_pert
	run_timestep_real(i_h_pert, i_u, i_v, o_h_pert, o_u, o_v, i_dt, i_simulation_timestamp);
}



bool SWE_Cart2D_TS_l_rexi::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	PDESWECart2DTS_BaseInterface::shackRegistration(io_shackDict);
	ts_l_direct.shackRegistration(io_shackDict);

	return error.exists();
}

void SWE_Cart2D_TS_l_rexi::run_timestep_real(
		const sweet::Data::Cart2D::DataSpectral &i_h_pert,	//!< prognostic variables
		const sweet::Data::Cart2D::DataSpectral &i_u,		//!< prognostic variables
		const sweet::Data::Cart2D::DataSpectral &i_v,		//!< prognostic variables

		sweet::Data::Cart2D::DataSpectral &o_h_pert,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &o_u,			//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &o_v,			//!< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
#if SWEET_BENCHMARK_TIMINGS
	sweet::Tools::StopwatchBox::getInstance().rexi_timestepping.start();
#endif

	final_timestep = false;

	if (rexi_use_direct_solution)
	{

		o_h_pert = i_h_pert;
		o_u = i_u;
		o_v = i_v;
		ts_l_direct.runTimestep(o_h_pert, o_u, o_v, i_dt, i_simulation_timestamp);

#if SWEET_BENCHMARK_TIMINGS
	sweet::Tools::StopwatchBox::getInstance().rexi_timestepping.stop();
#endif

		return;
	}

	if (i_dt <= 0)
		SWEETErrorFatal("Only constant time step size allowed");


	typedef std::complex<double> complex;

	std::size_t max_N = rexi_alphas.size();


	/*
	 * Request physical or spectral here to avoid parallel race conditions
	 */
/////#if !SWEET_USE_CART2D_SPECTRAL_SPACE
/////	i_h_pert.request_data_physical();
/////	i_u.request_data_physical();
/////	i_v.request_data_physical();
/////#else
/////	i_h_pert.request_data_spectral();
/////	i_u.request_data_spectral();
/////	i_v.request_data_spectral();
/////#endif

#if SWEET_MPI

#if SWEET_BENCHMARK_TIMINGS
	if (mpi_rank == 0)
		stopwatch_broadcast.start();
#endif

	std::size_t data_size = i_h_pert.cart2DDataConfig->spectral_array_data_number_of_elements;
	MPI_Bcast(i_h_pert.spectral_space_data, data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (std::isnan(i_h_pert.spectral_get(0,0).real()) || std::isnan(i_h_pert.spectral_get(0,0).imag()))
	{
		final_timestep = true;

#if SWEET_BENCHMARK_TIMINGS
	sweet::Tools::StopwatchBox::getInstance().rexi_timestepping.stop();
#endif
		return;
	}


	MPI_Bcast(i_u.spectral_space_data, data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(i_v.spectral_space_data, data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#if SWEET_BENCHMARK_TIMINGS
	if (mpi_rank == 0)
		stopwatch_broadcast.stop();
#endif

#endif



#if SWEET_THREADING_TIME_REXI
#	pragma omp parallel for schedule(static,1) default(none) shared(i_dt, i_h_pert, i_u, i_v, max_N, std::cout, std::cerr)
#endif
	for (int i = 0; i < num_local_rexi_par_threads; i++)
	{
#if SWEET_BENCHMARK_TIMINGS
		bool stopwatch_measure = false;
	#if SWEET_THREADING_TIME_REXI
		if (omp_get_thread_num() == 0)
	#endif
			if (mpi_rank == 0)
				stopwatch_measure = true;
#endif

#if SWEET_BENCHMARK_TIMINGS
		if (stopwatch_measure)
			stopwatch_preprocessing.start();
#endif

		double eta_bar = shackPDESWECart2D->h0;
		double g = shackPDESWECart2D->gravitation;

		sweet::Data::Cart2DComplex::Operators &opc = perThreadVars[i]->ops;

		sweet::Data::Cart2DComplex::DataSpectral &eta0 = perThreadVars[i]->eta0;
		sweet::Data::Cart2DComplex::DataSpectral &u0 = perThreadVars[i]->u0;
		sweet::Data::Cart2DComplex::DataSpectral &v0 = perThreadVars[i]->v0;

		sweet::Data::Cart2DComplex::DataSpectral &h_sum = perThreadVars[i]->h_sum;
		sweet::Data::Cart2DComplex::DataSpectral &u_sum = perThreadVars[i]->u_sum;
		sweet::Data::Cart2DComplex::DataSpectral &v_sum = perThreadVars[i]->v_sum;

		sweet::Data::Cart2DComplex::DataSpectral &eta = perThreadVars[i]->eta;


		h_sum.spectral_setZero();
		u_sum.spectral_setZero();
		v_sum.spectral_setZero();


#if !SWEET_USE_CART2D_SPECTRAL_SPACE
		eta0 = sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2DComplex_DataSpectral::grid_convert(i_h_pert);
		u0 = sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2DComplex_DataSpectral::grid_convert(i_u);
		v0 = sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2DComplex_DataSpectral::grid_convert(i_v);
#else

// TODO: find a nice solution for this
//		if (shackDict.rexi.use_half_poles)
		if (true)
		{
			eta0 = sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2DComplex_DataSpectral::grid_convert(i_h_pert);
			u0 = sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2DComplex_DataSpectral::grid_convert(i_u);
			v0 = sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2DComplex_DataSpectral::grid_convert(i_v);
		}
///		else
///		{
///			eta0 = sweet::Data::Cart2D::Convert::Cart2DDataSpectral_2_Cart2DDataSpectralComplex::spectral_convert(i_h_pert);
///			u0 = sweet::Data::Cart2D::Convert::Cart2DDataSpectral_2_Cart2DDataSpectralComplex::spectral_convert(i_u);
///			v0 = sweet::Data::Cart2D::Convert::Cart2DDataSpectral_2_Cart2DDataSpectralComplex::spectral_convert(i_v);
///		}
#endif

		/**
		 * SPECTRAL SOLVER - DO EVERYTHING IN SPECTRAL SPACE
		 *
		 * (alpha+dt L)U = U0
		 *
		 * (alpha/dt+L) (dt U) = U0
		 *
		 * (alpha/dt+L) U = U0/dt
		 */
		// convert to spectral space
		// scale with inverse of tau
		double inv_dt = (1.0/i_dt);
		eta0 = eta0*inv_dt;
		u0 = u0*inv_dt;
		v0 = v0*inv_dt;

#if SWEET_THREADING_TIME_REXI || SWEET_MPI

#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
		int local_thread_id = omp_get_thread_num();
#else
		int local_thread_id = 0;
#endif
		int global_thread_id = local_thread_id + num_local_rexi_par_threads*mpi_rank;

		std::size_t start = std::min(max_N, block_size*global_thread_id);
		std::size_t end = std::min(max_N, start+block_size);
#else
		std::size_t start = 0;
		std::size_t end = max_N;
#endif

		/*
		 * DO SUM IN PARALLEL
		 */
		//
		// precompute a bunch of values
		// this would belong to a serial part according to Amdahl's law
		//
		// (kappa + lhs_a)\eta = kappa/alpha*\eta_0 - (i_parameters.sim.f0*eta_bar/alpha) * rhs_b + rhs_a
		//
		sweet::Data::Cart2DComplex::DataSpectral rhs_a = eta_bar*(opc.diff_c_x(u0) + opc.diff_c_y(v0));
		sweet::Data::Cart2DComplex::DataSpectral rhs_b = (opc.diff_c_x(v0) - opc.diff_c_y(u0));

		sweet::Data::Cart2DComplex::DataSpectral lhs_a = (-g*eta_bar)*(perThreadVars[i]->ops.diff2_c_x + perThreadVars[i]->ops.diff2_c_y);

#if SWEET_BENCHMARK_TIMINGS
		if (stopwatch_measure)
			stopwatch_preprocessing.stop();
#endif

#if SWEET_BENCHMARK_TIMINGS
		if (stopwatch_measure)
			stopwatch_solve_rexi_terms.start();
#endif

		for (std::size_t n = start; n < end; n++)
		{
			// load alpha (a) and scale by inverse of tau
			complex alpha = -rexi_alphas[n]/i_dt;
			complex beta = rexi_betas[n];

			if (shackPDESWECart2D->cart2d_rotating_f0 == 0)
			{
				/*
				 * TODO: we can even get more performance out of this operations
				 * by partly using the real Fourier transformation
				 */
				sweet::Data::Cart2DComplex::DataSpectral rhs =
						eta0*alpha
						+ eta_bar*(opc.diff_c_x(u0) + opc.diff_c_y(v0))
					;

				sweet::Data::Cart2DComplex::DataSpectral lhs_a = (-g*eta_bar)*(perThreadVars[i]->ops.diff2_c_x + perThreadVars[i]->ops.diff2_c_y);
				sweet::Data::Cart2DComplex::DataSpectral lhs = lhs_a.spectral_addScalarAll(alpha*alpha);
				eta = rhs.spectral_div_element_wise(lhs);

				sweet::Data::Cart2DComplex::DataSpectral u1 = (u0 + g*opc.diff_c_x(eta))*(1.0/alpha);
				sweet::Data::Cart2DComplex::DataSpectral v1 = (v0 + g*opc.diff_c_y(eta))*(1.0/alpha);

				h_sum += eta*beta;
				u_sum += u1*beta;
				v_sum += v1*beta;
			}
			else
			{
				// load kappa (k)
				complex kappa = alpha*alpha + shackPDESWECart2D->cart2d_rotating_f0*shackPDESWECart2D->cart2d_rotating_f0;

				/*
				 * TODO: we can even get more performance out of this operations
				 * by partly using the real Fourier transformation
				 */
				sweet::Data::Cart2DComplex::DataSpectral rhs =
						(kappa/alpha) * eta0
						+ (-shackPDESWECart2D->cart2d_rotating_f0*eta_bar/alpha) * rhs_b
						+ rhs_a
					;

				sweet::Data::Cart2DComplex::DataSpectral lhs = lhs_a.spectral_addScalarAll(kappa);
				eta = rhs.spectral_div_element_wise(lhs);

				sweet::Data::Cart2DComplex::DataSpectral uh = u0 + g*opc.diff_c_x(eta);
				sweet::Data::Cart2DComplex::DataSpectral vh = v0 + g*opc.diff_c_y(eta);

				sweet::Data::Cart2DComplex::DataSpectral u1 = (alpha/kappa) * uh     - (shackPDESWECart2D->cart2d_rotating_f0/kappa) * vh;
				sweet::Data::Cart2DComplex::DataSpectral v1 = (shackPDESWECart2D->cart2d_rotating_f0/kappa) * uh + (alpha/kappa) * vh;

				sweet::Data::Cart2D::DataSpectral tmp(h_sum.cart2DDataConfig);

				h_sum += eta*beta;
				u_sum += u1*beta;
				v_sum += v1*beta;
			}
		}

#if SWEET_BENCHMARK_TIMINGS
		if (stopwatch_measure)
			stopwatch_solve_rexi_terms.stop();
#endif
	}




#if SWEET_BENCHMARK_TIMINGS
	if (mpi_rank == 0)
		stopwatch_reduce.start();
#endif

#if SWEET_THREADING_TIME_REXI


#if !SWEET_USE_CART2D_SPECTRAL_SPACE
	o_h_pert.grid_setZero();
	o_u.grid_setZero();
	o_v.grid_setZero();

	for (int n = 0; n < num_local_rexi_par_threads; n++)
	{
///		perThreadVars[n]->h_sum.request_data_physical();
///		perThreadVars[n]->u_sum.request_data_physical();
///		perThreadVars[n]->v_sum.request_data_physical();

		// sum real-valued elements
		#pragma omp parallel for schedule(static)
		for (std::size_t i = 0; i < i_h_pert.cart2DDataConfig->grid_number_elements; i++)
			o_h_pert.grid_space_data[i] += perThreadVars[n]->h_sum.grid_space_data[i].real();

		#pragma omp parallel for schedule(static)
		for (std::size_t i = 0; i < i_h_pert.cart2DDataConfig->grid_number_elements; i++)
			o_u.grid_space_data[i] += perThreadVars[n]->u_sum.grid_space_data[i].real();

		#pragma omp parallel for schedule(static)
		for (std::size_t i = 0; i < i_h_pert.cart2DDataConfig->grid_number_elements; i++)
			o_v.grid_space_data[i] += perThreadVars[n]->v_sum.grid_space_data[i].real();
	}
#else

	o_h_pert.spectral_setZero();
	o_u.spectral_setZero();
	o_v.spectral_setZero();

	for (int n = 0; n < num_local_rexi_par_threads; n++)
	{
///		perThreadVars[n]->h_sum.request_data_spectral();
///		perThreadVars[n]->u_sum.request_data_spectral();
///		perThreadVars[n]->v_sum.request_data_spectral();

		sweet::Data::Cart2D::DataSpectral tmp(perThreadVars[n]->ops.cart2DDataConfig);

		o_h_pert = o_h_pert + sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2D_DataSpectral::convert_real(perThreadVars[n]->h_sum);
		o_u = o_u + sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2D_DataSpectral::convert_real(perThreadVars[n]->u_sum);
		o_v = o_v + sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2D_DataSpectral::convert_real(perThreadVars[n]->v_sum);
	}
#endif


#else

#if !SWEET_USE_CART2D_SPECTRAL_SPACE
	o_h_pert = sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2D_DataSpectral::grid_convert(perThreadVars[0]->h_sum);
	o_u = sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2D_DataSpectral::grid_convert(perThreadVars[0]->u_sum);
	o_v = sweet::Data::Cart2D::Convert::DataSpectral_2_Cart2D_DataSpectral::grid_convert(perThreadVars[0]->v_sum);

#else

// TODO: find a nice solution for this
//		if (shackDict.rexi.use_half_poles)
	if (true)
	{
		o_h_pert = sweet::Data::Cart2DComplex::Convert::DataSpectral_2_Cart2D_DataSpectral::convert_real(perThreadVars[0]->h_sum);
		o_u = sweet::Data::Cart2DComplex::Convert::DataSpectral_2_Cart2D_DataSpectral::convert_real(perThreadVars[0]->u_sum);
		o_v = sweet::Data::Cart2DComplex::Convert::DataSpectral_2_Cart2D_DataSpectral::convert_real(perThreadVars[0]->v_sum);
	}
////	else
////	{
////		o_h_pert = sweet::Data::Cart2D::Convert::Cart2DDataSpectralComplex_2_Cart2DDataSpectral::spectral_convert_grid_real(perThreadVars[0]->h_sum);
////		o_u = sweet::Data::Cart2D::Convert::Cart2DDataSpectralComplex_2_Cart2DDataSpectral::spectral_convert_grid_real(perThreadVars[0]->u_sum);
////		o_v = sweet::Data::Convert::Cart2DDataSpectralComplex_2_Cart2DData::spectral_convert_grid_real(perThreadVars[0]->v_sum);
////	}
#endif


	if (rexi_gamma.real() != 0)
	{
		o_h_pert += rexi_gamma.real() * i_h_pert;
		o_u += rexi_gamma.real() * i_u;
		o_v += rexi_gamma.real() * i_v;
	}


#endif


#if SWEET_MPI

#if !SWEET_USE_CART2D_SPECTRAL_SPACE ||1
	sweet::Data::Cart2D::DataSpectral tmp(o_h_pert.cart2DDataConfig);

///	o_h_pert.request_data_physical();
	int retval = MPI_Reduce(o_h_pert.spectral_space_data, tmp.spectral_space_data, data_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (retval != MPI_SUCCESS)
	{
		std::cerr << "MPI FAILED!" << std::endl;
		exit(1);
	}

	std::swap(o_h_pert.spectral_space_data, tmp.spectral_space_data);

//	o_u.request_data_physical();
	MPI_Reduce(o_u.spectral_space_data, tmp.spectral_space_data, data_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	std::swap(o_u.spectral_space_data, tmp.spectral_space_data);

///	o_v.request_data_physical();
	MPI_Reduce(o_v.spectral_space_data, tmp.spectral_space_data, data_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	std::swap(o_v.spectral_space_data, tmp.spectral_space_data);

#else

	#error "TODO: spectral version"

#endif

#endif


#if SWEET_BENCHMARK_TIMINGS
	if (mpi_rank == 0)
		stopwatch_reduce.stop();
#endif

#if SWEET_BENCHMARK_TIMINGS
	sweet::Tools::StopwatchBox::getInstance().rexi_timestepping.stop();
#endif
}



void SWE_Cart2D_TS_l_rexi::runTimestep(
		sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	runTimestep(
			io_h, io_u, io_v,
			io_h, io_u, io_v,
			i_dt,
			i_simulation_timestamp
		);
}


void SWE_Cart2D_TS_l_rexi::cleanup()
{

	for (std::vector<PerThreadVars*>::iterator iter = perThreadVars.begin(); iter != perThreadVars.end(); ++iter)
	{
		PerThreadVars* p = *iter;
		delete p;
	}

	perThreadVars.resize(0);
}



void SWE_Cart2D_TS_l_rexi::MPI_quitWorkers(
		sweet::Data::Cart2D::Config *i_cart2DDataConfig
)
{
#if SWEET_MPI
	sweet::Data::Cart2D::DataSpectral dummyData(i_cart2DDataConfig);
	dummyData.spectral_setValue(NAN);

	MPI_Bcast(dummyData.spectral_space_data, dummyData.cart2DDataConfig->spectral_array_data_number_of_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}




SWE_Cart2D_TS_l_rexi::~SWE_Cart2D_TS_l_rexi()
{
	cleanup();
}

