/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONSPHERE2D_BENCHMARKS_SHACKPDEADVECTIONSPHERE2DBENCHMARKS_HPP
#define PROGRAMS_PDE_ADVECTIONSPHERE2D_BENCHMARKS_SHACKPDEADVECTIONSPHERE2DBENCHMARKS_HPP

#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Vector/Vector.hpp>
#include <sweet/Shacks/Base.hpp>
#include <string>
#include <iostream>

/**
 * Values and parameters to setup benchmarks simulations
 */
class ShackPDEAdvectionSphere2DBenchmarks	:
		public sweet::Shacks::Base
{
public:
	//! seed for random number generator
	int random_seed = 0;

	//! benchmark scenario
	std::string benchmark_name = "";

	//! rotation angle for advection equation
	double sphere2d_advection_rotation_angle = 0;


	/*
	 * Get updated velocities for particular point in time
	 */
	void (*getVelocities)(
			sweet::Data::Sphere2D::DataGrid&,
			sweet::Data::Sphere2D::DataGrid&,
			double i_time,
			ShackPDEAdvectionSphere2DBenchmarks* user_ptr
	) = nullptr;
	ShackPDEAdvectionSphere2DBenchmarks *getVelocitiesUserData = nullptr;

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


	void printProgramArguments(const std::string& i_prefix = "") override
	{
		std::cout << i_prefix << std::endl;
		std::cout << i_prefix << "SIMULATION SETUP PARAMETERS:" << std::endl;
		std::cout << i_prefix << "	--random-seed [int]		random seed for random number generator" << std::endl;
		std::cout << i_prefix << "	--benchmark-name [string]	benchmark name" << std::endl;
		std::cout << i_prefix << "	--advection-velocity=[float],[float],[float]	advection velocity components (x, y, rotational)" << std::endl;
		std::cout << i_prefix << "	--benchmark-override-simvars [bool]	Allow overwriting simulation variables by benchmark (default: 1)" << std::endl;
		std::cout << i_prefix << "	--benchmark-setup-dealiased [bool]	Use dealiasing for setup (default: 1)" << std::endl;
		std::cout << i_prefix << "	-x [float]				x coordinate for setup \\in [0;1], default=0.5" << std::endl;
		std::cout << i_prefix << "	-y [float]				y coordinate for setup \\in [0;1], default=0.5" << std::endl;
		std::cout << i_prefix << "	-r [radius]				scale factor of radius for initial condition, default=1" << std::endl;
		std::cout << i_prefix << "	--initial-coord-x [float]		Same as -x" << std::endl;
		std::cout << i_prefix << "	--initial-coord-y [float]		Same as -y" << std::endl;
		std::cout << i_prefix << "	--benchmark-advection-rotation-angle [float]	Rotation angle for e.g. advection test case" << std::endl;
		std::cout << i_prefix << std::endl;
	}

	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa) override
	{
		i_pa.getArgumentValueByKey("--random-seed", random_seed);
		i_pa.getArgumentValueByKey("--benchmark-name", benchmark_name);

		i_pa.getArgumentValueByKey("--benchmark-advection-rotation-angle", sphere2d_advection_rotation_angle);

		if (error.exists())
			return error.forwardWithPositiveReturn(i_pa.error);

		if (random_seed >= 0)
			srandom(random_seed);

		return error.forwardWithPositiveReturn(i_pa.error);
	}

	void printShack(
		const std::string& i_prefix = ""
	) override
	{
		std::cout << i_prefix << std::endl;
		std::cout << i_prefix << "BENCHMARK:" << std::endl;
		std::cout << i_prefix << " + random_seed: " << random_seed << std::endl;
		std::cout << i_prefix << " + benchmark_name: " << benchmark_name << std::endl;
		std::cout << i_prefix << " + sphere2d_advection_rotation_angle: " << sphere2d_advection_rotation_angle << std::endl;
		std::cout << i_prefix << " + input_data_filenames:" << std::endl;
		std::cout << i_prefix << std::endl;
	}
};


#endif

