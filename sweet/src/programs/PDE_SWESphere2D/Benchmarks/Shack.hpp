/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_SHACK_HPP
#define PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_SHACK_HPP

#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Vector/Vector.hpp>
#include <sweet/Shacks/Base.hpp>
#include <string>
#include <iostream>

#include "../DataContainer/Topography.hpp"

namespace PDE_SWESphere2D {
namespace Benchmarks {

/**
 * Values and parameters to setup benchmarks simulations
 */
class Shack	:
		public sweet::Shacks::Base
{
public:
	//! seed for random number generator
	int random_seed = 0;

	//! benchmark scenario
	std::string benchmark_name = "";

	//! rotation angle for advection equation
	double benchmark_sphere2d_advection_rotation_angle = 0;

	/**
	 * Flag to indicate the presence of topography
	 */
	bool use_topography = false;

	/*
	 * Topography itself
	 */
	//////sweet::Data::Sphere2D::DataGrid h_topography;
	PDE_SWESphere2D::DataContainer::Topography topography;


	//! Galewsky-benchmark specific: velocity
	double benchmark_galewsky_umax = -1;

	//! Galewsky-benchmark specific: amplitude of bump
	double benchmark_galewsky_hamp = -1;

	//! Galewsky-benchmark specific: latitude coordinate
	double benchmark_galewsky_phi2 = -1;

	std::string benchmark_galewsky_geostrophic_setup = "analytical";


	/*
	 * Get updated velocities for particular point in time
	 */
	void (*getVelocities)(
			sweet::Data::Sphere2D::DataGrid&,
			sweet::Data::Sphere2D::DataGrid&,
			double i_time,
			Shack* user_ptr
	) = nullptr;
	Shack *getVelocitiesUserData = nullptr;

	/**
	 * Callback for special benchmark
	 */
	void (*callback_slComputeDeparture3rdOrder)(
			void *i_this,
			const sweet::Data::Vector::Vector<double> &i_pos_lon_A,	//!< longitude coordinate to compute the velocity for
			const sweet::Data::Vector::Vector<double> &i_pos_lat_A,	//!< latitude coordinate to compute the velocity for
			sweet::Data::Vector::Vector<double> &o_pos_lon_D,		//!< velocity along longitude
			sweet::Data::Vector::Vector<double> &o_pos_lat_D,		//!< velocity along latitude
			double i_dt,
			double i_timestamp_arrival			//!< timestamp at arrival point
	);
	void *slComputeDeparture3rdOrderUserData = nullptr;



	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << std::endl;
		std::cout << i_prefix << "SIMULATION SETUP PARAMETERS:" << std::endl;
		std::cout << i_prefix << "	--benchmark-random-seed [int]		random seed for random number generator" << std::endl;
		std::cout << i_prefix << "	--benchmark-name [string]	benchmark name" << std::endl;
		std::cout << i_prefix << "	--benchmark-override-simvars [bool]	Allow overwriting simulation variables by benchmark (default: 1)" << std::endl;
		std::cout << i_prefix << "	--benchmark-setup-dealiased [bool]	Use dealiasing for setup (default: 1)" << std::endl;
		std::cout << i_prefix << "	--benchmark-advection-rotation-angle [float]	Rotation angle for e.g. advection test case" << std::endl;
		std::cout << i_prefix << "	--benchmark-galewsky-geostropic-setup [str]	'analytical'/'numerical' setup of geostropic balance" << std::endl;
		std::cout << i_prefix << std::endl;
	}

	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--benchmark-random-seed", random_seed);
		i_pa.getArgumentValueByKey("--benchmark-name", benchmark_name);

		i_pa.getArgumentValueByKey("--benchmark-advection-rotation-angle", benchmark_sphere2d_advection_rotation_angle);

		i_pa.getArgumentValueByKey("--benchmark-galewsky-umax", benchmark_galewsky_umax);
		i_pa.getArgumentValueByKey("--benchmark-galewsky-hamp", benchmark_galewsky_hamp);
		i_pa.getArgumentValueByKey("--benchmark-galewsky-phi2", benchmark_galewsky_phi2);
		i_pa.getArgumentValueByKey("--benchmark-galewsky-geostropic-setup", benchmark_galewsky_geostrophic_setup);

		if (error.exists())
			return error.forwardWithPositiveReturn(i_pa.error);

		if (random_seed >= 0)
			srandom(random_seed);

		return error.forwardWithPositiveReturn(i_pa.error);
	}

	void printShack(
		const std::string& i_prefix = ""
	)
	{
		std::cout << i_prefix << std::endl;
		std::cout << i_prefix << "BENCHMARK:" << std::endl;
		std::cout << i_prefix << " + benchmark_name: " << benchmark_name << std::endl;
		std::cout << i_prefix << " + benchmark_random_seed: " << random_seed << std::endl;
		std::cout << i_prefix << " + benchmark_sphere2d_advection_rotation_angle: " << benchmark_sphere2d_advection_rotation_angle << std::endl;
		std::cout << i_prefix << " + benchmark_galewsky_umax: " << benchmark_galewsky_umax << std::endl;
		std::cout << i_prefix << " + benchmark_galewsky_hamp: " << benchmark_galewsky_hamp << std::endl;
		std::cout << i_prefix << " + benchmark_galewsky_phi2: " << benchmark_galewsky_phi2 << std::endl;
		std::cout << i_prefix << " + benchmark_galewsky_geostrophic_setup: " << benchmark_galewsky_geostrophic_setup << std::endl;
		std::cout << i_prefix << std::endl;
	}
};

}}

#endif
