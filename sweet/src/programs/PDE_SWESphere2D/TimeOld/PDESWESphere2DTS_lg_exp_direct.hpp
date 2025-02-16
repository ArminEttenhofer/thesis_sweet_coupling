/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LG_EXP_DIRECT_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LG_EXP_DIRECT_HPP



#include <complex>
#include <string.h>
#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Data/Sphere2DComplex/DataSpectral.hpp>
#include <sweet/Data/Sphere2DComplex/Operators.hpp>
#include <sweet/ExpIntegration/ExpFunction.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include "../TimeHelpers/SWERexiTerm_SPH.hpp"

#include "PDESWESphere2DTS_BaseInterface.hpp"


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



class PDESWESphere2DTS_lg_exp_direct	: public PDESWESphere2DTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

public:
	bool setup_main(
			const sweet::Data::Sphere2D::Operators *io_ops,
			const std::string &i_function_name
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;

	std::string getIDString() override;


private:
	typedef std::complex<double> complex;

private:

	const sweet::Data::Sphere2D::Config *sphere2DDataConfig;

public:
	sweet::ExpIntegration::ExpFunction<double> expFunctions;


private:
	/*
	 * Function name to be used by REXI
	 */
	std::string function_name;


private:
	void reset();


public:

	void run_timestep_lg_exp(
		sweet::Data::Sphere2D::DataSpectral &io_prog_phi,
		sweet::Data::Sphere2D::DataSpectral &io_prog_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_prog_div,

		double i_dt,
		double i_simulation_timestamp
	);

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

	PDESWESphere2DTS_lg_exp_direct();

	virtual ~PDESWESphere2DTS_lg_exp_direct();
};


#endif
