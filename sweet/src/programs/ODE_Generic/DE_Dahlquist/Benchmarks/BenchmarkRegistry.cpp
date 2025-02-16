/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "BenchmarkRegistry.hpp"
#include "Benchmark_InitToOne.hpp"

namespace ODE_Generic {
namespace DE_Dahlquist {
namespace Benchmarks {

BenchmarkRegistry::BenchmarkRegistry()	:
	shackDict(nullptr),
	benchmark(nullptr)
{
}


bool BenchmarkRegistry::setup_1_registerAllBenchmark()
{
	_registered_benchmarks.push_back(static_cast<Base*>(new Benchmark_InitToOne));

	return true;
}



bool BenchmarkRegistry::setup_2_shackRegistration(
	sweet::Shacks::Dictionary *io_shackDict
)
{
	shackDict = io_shackDict;

	shackBenchmarks = shackDict->getAutoRegistration<Benchmarks::Shack>();

	for (std::size_t i = 0; i < _registered_benchmarks.size(); i++)
	{
		_registered_benchmarks[i]->shackRegistration(io_shackDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_registered_benchmarks[i]);
	}

	return true;
}



bool BenchmarkRegistry::setup_3_benchmarkDetection(
		const std::string &i_benchmark_name
)
{
	SWEET_ASSERT(benchmark == nullptr);

	std::string benchmark_name;

	if (i_benchmark_name != "")
	{
		benchmark_name = i_benchmark_name;
	}
	else
	{
		if (shackBenchmarks->benchmark_name == "")
		{
			printAvailableBenchmarks();
			return error.set("Please choose benchmark with --benchmark-name=...");
		}

		benchmark_name = shackBenchmarks->benchmark_name;
	}


	/*
	 * Find right one
	 */
	for (std::size_t i = 0; i < _registered_benchmarks.size(); i++)
	{
		Base *ts = _registered_benchmarks[i];

		if (ts->implements_benchmark(benchmark_name))
		{
			if (benchmark != nullptr)
			{
				//std::cout << "Processing " << i+1 << "th element" << std::endl;
				error.set("Duplicate implementation for benchmark "+benchmark_name);
			}

			//std::cout << "Benchmark detection: found match with benchmark id " << i+1 << std::endl;
			benchmark = ts;
		}
	}

	if (benchmark == nullptr)
	{
		printAvailableBenchmarks();
		return error.set("No valid --benchmark-name '"+benchmark_name+"' provided");
	}

	return true;
}


bool BenchmarkRegistry::setup_4_benchmarkSetup()
{
	SWEET_ASSERT(benchmark != nullptr);

	benchmark->setup_1_shackData();
	ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(*benchmark);

	return true;
}


void BenchmarkRegistry::clear_3_benchmarkDetection()
{
	benchmark = nullptr;
}


void BenchmarkRegistry::clear()
{
	clear_3_benchmarkDetection();

	shackDict = nullptr;

	_benchmarksFreeAll();
}



void BenchmarkRegistry::_benchmarksFreeAll(
		Base *skip_this
)
{
	for (std::size_t i = 0; i < _registered_benchmarks.size(); i++)
	{
		Base *ts = _registered_benchmarks[i];

		if (ts == skip_this)
			continue;

		delete ts;
	}

	_registered_benchmarks.clear();
	SWEET_ASSERT(_registered_benchmarks.size() == 0);
}


void BenchmarkRegistry::printAvailableBenchmarks(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << "Benchmarks (START)" << std::endl;
	o_ostream << "********************************************************************************" << std::endl;

	for (std::size_t i = 0; i < _registered_benchmarks.size(); i++)
	{
		std::cout << _registered_benchmarks[i]->printHelp();
		std::cout << std::endl;
	}

	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << "Benchmarks (END)" << std::endl;
	o_ostream << "********************************************************************************" << std::endl;

}


BenchmarkRegistry::~BenchmarkRegistry()
{
	clear();
}

}}}
