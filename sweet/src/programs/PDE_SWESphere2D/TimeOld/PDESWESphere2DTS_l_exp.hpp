/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_L_EXP_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_L_EXP_HPP


#ifndef SWEET_MPI
#define SWEET_MPI 1
#endif


#include <complex>
#include <sweet/ExpIntegration/REXI/REXI_Terry.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <string.h>
#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2DComplex/DataSpectral.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Data/Sphere2DComplex/Operators.hpp>
#include <sweet/ExpIntegration/ExpFunction.hpp>
#include "../TimeHelpers/SWERexiTerm_SPH.hpp"
#include "PDESWESphere2DTS_BaseInterface.hpp"
#include "PDESWESphere2DTS_l_exp_direct_special.hpp"
#include "PDESWESphere2DTS_lg_exp_lc_taylor.hpp"


#ifndef SWEET_BENCHMARK_TIMINGS
	#define SWEET_BENCHMARK_TIMINGS 1
#endif

#ifndef SWEET_REXI_TIMINGS_ADDITIONAL_BARRIERS
	#define SWEET_REXI_TIMINGS_ADDITIONAL_BARRIERS 1
#endif



#if SWEET_BENCHMARK_TIMINGS
#include <sweet/Tools/Stopwatch.hpp>
#endif

#if SWEET_MPI
	#include <mpi.h>
#endif



class PDESWESphere2DTS_l_exp	: public PDESWESphere2DTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

public:
	bool setup_variant_10(
			const sweet::Data::Sphere2D::Operators *io_ops,
			sweet::ExpIntegration::Shack *i_shackExpIntegration,
			const std::string &i_function_name,
			double i_timestepSize,
			bool i_use_f_sphere2D,
			bool i_no_coriolis
	);

public:
	bool setup_variant_50(
			const sweet::Data::Sphere2D::Operators *io_ops,
			sweet::ExpIntegration::Shack *i_shackExpIntegration,
			const std::string &i_function_name,
			double i_timestepSize,
			bool i_use_f_sphere2D,
			bool i_no_coriolis,
			int i_timestepping_order
	);

public:
	bool setup_variant_100(
			const sweet::Data::Sphere2D::Operators *io_ops,
			sweet::ExpIntegration::Shack *i_shackExpIntegration,
			const std::string &i_function_name,
			const std::string &i_exp_method,
			double i_timestepSize,
			bool i_use_f_sphere2D,
			bool i_no_coriolis,
			int i_timestepping_order,
			bool i_use_rexi_sphere2d_solver_preallocation
	);


public:
	bool shackRegistration(sweet::Shacks::Dictionary *io_shackDict) override;

	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;
	std::string getIDString() override;


private:
	typedef std::complex<double> complex;


public:
	std::vector<std::complex<double>> _rexi_alphas;
	std::vector<std::complex<double>> _rexi_betas;
	std::complex<double> _rexi_gamma;

	sweet::ExpIntegration::ExpFunction<double> _expFunction;


private:
	/*!
	 * Time step size of REXI
	 */
	double timestep_size;

	/*!
	 * Function name to be used by REXI
	 */
	std::string function_name;

	/*!
	 * Exponential integration method to use
	 */
	std::string exp_method;

	/*!
	 * Don't use any Coriolis effect (reduction to very simple Helmholtz problem)
	 */
	bool no_coriolis;

	/*!
	 * Assume f-sphere2D (reduction to Helmholtz problem)
	 */
	bool use_f_sphere2D;

	/*!
	 * Preallocate the REXI matrices
	 */
	bool use_rexi_sphere2d_solver_preallocation;

	/*!
	 * True, if EXP method is "direct"
	 */
	bool use_exp_method_direct_solution;

	/*!
	 * True, if EXP method is "ss_taylor"
	 */
	bool use_exp_method_strang_split_taylor;

	/*!
	 * True, if a REXI method is used
	 */
	bool use_exp_method_rexi;


	std::size_t block_size;

	class PerThreadVars
	{
	public:
		std::vector<SWERexiTerm_SPH> sweREXISolvers;

		std::vector< std::complex<double> > alpha;
		std::vector< std::complex<double> > beta;

		sweet::Data::Sphere2D::DataSpectral accum_phi;
		sweet::Data::Sphere2D::DataSpectral accum_vrt;
		sweet::Data::Sphere2D::DataSpectral accum_div;
	};

	// per-thread allocated variables to avoid NUMA domain effects
	std::vector<PerThreadVars*> perThreadVars;

	// number of threads to be used
	int num_local_rexi_par_threads;

	// number of threads to be used
	int num_global_threads;


#if SWEET_MPI
	// MPI communicator
	MPI_Comm mpi_comm;

	// number of mpi ranks to be used
	int mpi_rank;

	// MPI ranks
	int num_mpi_ranks;
#endif

	PDESWESphere2DTS_l_exp_direct_special timestepping_method_l_exp_direct_special;

	PDESWESphere2DTS_lg_exp_lc_taylor timestepping_method_lg_exp_lc_exp;


private:
	void reset();


public:
	PDESWESphere2DTS_l_exp();

private:
	void p_update_coefficients();

	void p_get_workload_start_end(
			std::size_t &o_start,
			std::size_t &o_end,
			int i_local_thread_id
	);

public:
	void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vort,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;


	void runTimestep(
			const sweet::Data::Sphere2D::DataSpectral &i_h,
			const sweet::Data::Sphere2D::DataSpectral &i_u,
			const sweet::Data::Sphere2D::DataSpectral &i_v,

			sweet::Data::Sphere2D::DataSpectral &o_h,
			sweet::Data::Sphere2D::DataSpectral &o_u,
			sweet::Data::Sphere2D::DataSpectral &o_v,

			double i_fixed_dt,
			double i_simulation_timestamp
	);


	/**
	 * Solve the REXI of \f$ U(t) = exp(L*t) \f$
	 *
	 * See
	 * 		doc/rexi/understanding_rexi.pdf
	 * for further information
	 */
public:
	bool run_timestep_rexi_vectorinvariant_progphivortdiv(
		sweet::Data::Sphere2D::DataSpectral &io_phi0,
		sweet::Data::Sphere2D::DataSpectral &io_u0,
		sweet::Data::Sphere2D::DataSpectral &io_v0,

		double i_timestepSize,	//!< timestep size

		const sweet::Shacks::Dictionary &i_parameters
	);


	virtual ~PDESWESphere2DTS_l_exp();
};


#endif
