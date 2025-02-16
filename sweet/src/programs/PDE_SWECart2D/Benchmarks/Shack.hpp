/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWECART2D_BENCHMARKS_SHACKPDESWECART2DBENCHMARKS_HPP
#define PROGRAMS_PDE_SWECART2D_BENCHMARKS_SHACKPDESWECART2DBENCHMARKS_HPP

#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Shacks/Base.hpp>
#include <string>
#include <iostream>
#include <sweet/Tools/ProgramArguments.hpp>



/**
 * simulation coefficients
 */

namespace PDE_SWECart2D {
namespace Benchmarks {

class Shack	:
		public sweet::Shacks::Base
{
public:

	//! seed for random number generator
	int random_seed = 0;

	//! benchmark scenario
	std::string benchmark_name = "";

	//! Normal modes benchmark scenario
	std::string benchmark_normal_modes_case = "";

	//! radius
	double object_scale = 1;

	//! setup coordinate of e.g. radial breaking dam, x-placement \in [0;1]
	double object_coord_x = 0.5;

	//! setup coordinate of e.g. radial breaking dam, y-placement \in [0;1]
	double object_coord_y = 0.5;

	//! load external forces if available from benchmark scenario
	void (*getExternalForcesCallback)(int, double, sweet::Data::Cart2D::DataSpectral*, PDE_SWECart2D::Benchmarks::Shack*) = nullptr;// = &fun_no_forces;		//! SET TO NULLPTR
	PDE_SWECart2D::Benchmarks::Shack *getExternalForcesUserData = nullptr;


	/**
	 * Velocity and additional parameter for advection test cases
	 */
	double advection_velocity[3] = {0, 0};


	void printProgramArguments(const std::string& i_prefix = "") override
	{
		std::cout << i_prefix << std::endl;
		std::cout << i_prefix << "SIMULATION SETUP PARAMETERS:" << std::endl;
		std::cout << i_prefix << "	--random-seed [int]		random seed for random number generator" << std::endl;
		std::cout << i_prefix << "	--benchmark-name [string]	benchmark name" << std::endl;
		std::cout << i_prefix << "	-x [float]				x coordinate for setup \\in [0;1], default=0.5" << std::endl;
		std::cout << i_prefix << "	-y [float]				y coordinate for setup \\in [0;1], default=0.5" << std::endl;
		std::cout << i_prefix << "	-r [radius]				scale factor of radius for initial condition, default=1" << std::endl;
		std::cout << i_prefix << "	--advection-velocity=[float],[float],[float]	advection velocity components (x, y, rotational)" << std::endl;

		std::cout << "" << std::endl;
	}

	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa) override
	{
		i_pa.getArgumentValueByKey("--random-seed", random_seed);
		i_pa.getArgumentValueBy2Keys("--initial-coord-x", "-x", object_coord_x);
		i_pa.getArgumentValueBy2Keys("--initial-coord-y", "-y", object_coord_y);
		i_pa.getArgumentValueByKey("--benchmark-name", benchmark_name);

		i_pa.getArgumentValueByKey("--benchmark-normal-modes-case", benchmark_normal_modes_case);
		i_pa.getArgumentValueByKey("-r", object_scale);

		std::string tmp;
		if (i_pa.getArgumentValueByKey("--advection-velocity", tmp))
			sweet::Tools::StringSplit::split3double(tmp, &advection_velocity[0], &advection_velocity[1], &advection_velocity[2]);


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
		std::cout << std::endl;
		std::cout << "BENCHMARK:" << std::endl;
		std::cout << " + random_seed: " << random_seed << std::endl;
		std::cout << " + benchmark_name: " << benchmark_name << std::endl;
		std::cout << " + benchmark_normal_modes_case: " << benchmark_normal_modes_case << std::endl;
		std::cout << " + object_scale: " << object_scale << std::endl;
		std::cout << " + object_coord_x: " << object_coord_x << std::endl;
		std::cout << " + object_coord_y: " << object_coord_y << std::endl;
		std::cout << i_prefix << " + advection_velocity (x, y, rotation speed): " << advection_velocity[0] << ", " << advection_velocity[1] << std::endl;
		std::cout << std::endl;
	}
};


}}



#endif
