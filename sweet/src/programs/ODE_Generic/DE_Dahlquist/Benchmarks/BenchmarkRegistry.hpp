#ifndef PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_BENCHMARKS_BENCHMARKREGISTRY_HPP
#define PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_BENCHMARKS_BENCHMARKREGISTRY_HPP


#include <iostream>

#include <sweet/Error/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>

#include "Base.hpp"
#include "Shack.hpp"

namespace ODE_Generic {
namespace DE_Dahlquist {
namespace Benchmarks {


class BenchmarkRegistry
{
	std::vector<Base*> _registered_benchmarks;

	sweet::Shacks::Dictionary *shackDict;
	Shack *shackBenchmarks;

public:
	sweet::Error::Base error;

	Base *benchmark = nullptr;

	BenchmarkRegistry();

	bool setup_1_registerAllBenchmark();

	bool setup_2_shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
	);

	bool setup_3_benchmarkDetection(
			const std::string &i_benchmark_name = ""
	);

	bool setup_4_benchmarkSetup();

	void clear_3_benchmarkDetection();

	void clear();

private:
	void _benchmarksFreeAll(
			Base *skip_this = nullptr
	);

public:
	void printAvailableBenchmarks(
		std::ostream &o_ostream = std::cout,
		const std::string &i_prefix = ""
	);

public:
	~BenchmarkRegistry();
};

}}}

#endif
