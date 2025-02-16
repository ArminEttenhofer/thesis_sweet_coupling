/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef PROGRAMS_ODE_SCALAR_ODESCALARBENCHMARKSCOMBINED_HPP
#define PROGRAMS_ODE_SCALAR_ODESCALARBENCHMARKSCOMBINED_HPP

#include <programs/ODE_Scalar/benchmarks/ODEScalarBench_Default.hpp>
#include <programs/ODE_Scalar/benchmarks/ShackODEScalarBenchmarks.hpp>
#include <iostream>
#include <sweet/Error/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>


class ODEScalarBenchmarksCombined
{
public:
	sweet::Error::Base error;
	ShackODEScalarBenchmarks *shackBenchmarks;

	ODEScalarBenchmarksCombined()	:
		shackBenchmarks(nullptr)
	{

	}


	/*
	 * Special function to register shacks for benchmarks.
	 *
	 * This is in particular important for the --help output function to include all benchmarks.
	 */
	bool shackRegistration(
			sweet::Shacks::Dictionary &io_shackDict
	)
	{
		shackBenchmarks = io_shackDict.getAutoRegistration<ShackODEScalarBenchmarks>();

		return error.forwardWithPositiveReturn(io_shackDict.error);
	}

	void clear()
	{
		shackBenchmarks = nullptr;
	}


public:
	template <typename T>
	bool setupInitialConditions(
			T &o_u,
			sweet::Shacks::Dictionary &io_shackDict
	)
	{
		if (!shackBenchmarks->validateNonStationaryODE())
			return error.forwardWithPositiveReturn(shackBenchmarks->error);

		if (shackBenchmarks->benchmark_name == "")
			return error.set("ODEScalarBenchmarksCombined: Benchmark name not given, use --benchmark-name=[name]");

		if (shackBenchmarks->benchmark_name == "default")
		{
			ODEScalarBenchDefault default_bench;

			default_bench.shackRegistration(&io_shackDict);

			default_bench.setup();

			default_bench.setupBenchmark(
					o_u
			);

			return true;
		}

		shackBenchmarks->printProgramArguments();

		return error.set(std::string("Benchmark '")+shackBenchmarks->benchmark_name+ "' not found (or not available)");
	}
};



#endif
