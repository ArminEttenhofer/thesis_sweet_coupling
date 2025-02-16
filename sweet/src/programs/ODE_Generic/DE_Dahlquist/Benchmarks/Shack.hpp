/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_BENCHMARKS_SHACK_HPP
#define PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_BENCHMARKS_SHACK_HPP

#include <sweet/Shacks/Base.hpp>
#include <string>
#include <iostream>

namespace ODE_Generic {
namespace DE_Dahlquist {
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


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << std::endl;
		std::cout << i_prefix << "SIMULATION SETUP PARAMETERS:" << std::endl;
		std::cout << i_prefix << "	--benchmark-random-seed [int]		random seed for random number generator" << std::endl;
		std::cout << i_prefix << "	--benchmark-name [string]	benchmark name" << std::endl;
		std::cout << i_prefix << std::endl;
	}

	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--benchmark-random-seed", random_seed);
		i_pa.getArgumentValueByKey("--benchmark-name", benchmark_name);

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
		std::cout << i_prefix << " + benchmark_random_seed: " << random_seed << std::endl;
		std::cout << i_prefix << " + benchmark_name: " << benchmark_name << std::endl;
		std::cout << i_prefix << std::endl;
	}
};

}}}

#endif

