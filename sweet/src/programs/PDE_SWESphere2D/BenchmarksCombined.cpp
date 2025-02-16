/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "BenchmarksCombined.hpp"
#include "Benchmarks/zero.hpp"
#include "Benchmarks/galewsky.hpp"
#include "Benchmarks/gaussian_bump.hpp"
#include "Benchmarks/gaussian_bumps_pvd.hpp"
#include "Benchmarks/gaussian_bumps_test_cases.hpp"
#include "Benchmarks/three_gaussian_bumps.hpp"
#include "Benchmarks/williamson_1_advection_cos_bell.hpp"
#include "Benchmarks/williamson_1_advection_gauss_bump.hpp"
#include "Benchmarks/williamson_2_geostrophic_balance.hpp"
#include "Benchmarks/williamson_2_geostrophic_balance_linear.hpp"
#include "Benchmarks/williamson_2_geostrophic_balance_topography.hpp"
#include "Benchmarks/williamson_5_flow_over_mountain.hpp"
#include "Benchmarks/williamson_6_rossby_haurwitz_wave.hpp"
#include "Benchmarks/barotropic_vort_modes.hpp"
#include "Benchmarks/column.hpp"


namespace PDE_SWESphere2D {
namespace Benchmarks {

BenchmarksCombined::BenchmarksCombined()	:
	shackDict(nullptr),
	ops(nullptr),
	benchmark(nullptr)
{
}


bool BenchmarksCombined::setup_1_registerAllBenchmark()
{
	_registered_benchmarks.push_back(static_cast<BaseInterface*>(new gaussian_bumps_test_cases));
	_registered_benchmarks.push_back(static_cast<BaseInterface*>(new gaussian_bump));
	_registered_benchmarks.push_back(static_cast<BaseInterface*>(new gaussian_bumps_pvd));
	_registered_benchmarks.push_back(static_cast<BaseInterface*>(new galewsky));
	_registered_benchmarks.push_back(static_cast<BaseInterface*>(new three_gaussian_bumps));
	_registered_benchmarks.push_back(static_cast<BaseInterface*>(new williamson_1_advection_cos_bell));
	_registered_benchmarks.push_back(static_cast<BaseInterface*>(new williamson_1_advection_gauss_bump));
	_registered_benchmarks.push_back(static_cast<BaseInterface*>(new williamson_2_geostrophic_balance));
	_registered_benchmarks.push_back(static_cast<BaseInterface*>(new williamson_2_geostrophic_balance_linear));
	_registered_benchmarks.push_back(static_cast<BaseInterface*>(new williamson_2_geostrophic_balance_topography));
	_registered_benchmarks.push_back(static_cast<BaseInterface*>(new williamson_5_flow_over_mountain));
	_registered_benchmarks.push_back(static_cast<BaseInterface*>(new williamson_6_rossby_haurwitz_wave));
	_registered_benchmarks.push_back(static_cast<BaseInterface*>(new barotropic_vort_modes));
    _registered_benchmarks.push_back(static_cast<BaseInterface*>(new column));
    _registered_benchmarks.push_back(static_cast<BaseInterface*>(new zero));

	return true;
}



bool BenchmarksCombined::setup_2_shackRegistration(
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



bool BenchmarksCombined::setup_3_benchmarkDetection(
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
		BaseInterface *ts = _registered_benchmarks[i];

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


bool BenchmarksCombined::setup_4_benchmarkSetup_1_withoutOps()
{
	SWEET_ASSERT(benchmark != nullptr);

	benchmark->setup_1_shackData();
	ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(*benchmark);

	return true;
}


bool BenchmarksCombined::setup_5_benchmarkSetup_2_withOps(
		sweet::Data::Sphere2D::Operators* io_ops
)
{
	ops = io_ops;

	benchmark->setup_2_withOps(ops);
	ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(*benchmark);
}

void BenchmarksCombined::clear_3_benchmarkDetection()
{
	benchmark = nullptr;
}


void BenchmarksCombined::clear()
{
	clear_3_benchmarkDetection();

	shackDict = nullptr;
	ops = nullptr;

	_benchmarksFreeAll();
}



void BenchmarksCombined::_benchmarksFreeAll(
		BaseInterface *skip_this
)
{
	for (std::size_t i = 0; i < _registered_benchmarks.size(); i++)
	{
		BaseInterface *ts = _registered_benchmarks[i];

		if (ts == skip_this)
			continue;

		delete ts;
	}

	_registered_benchmarks.clear();
}


void BenchmarksCombined::printAvailableBenchmarks(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << "Benchmarks (START)" << std::endl;
	o_ostream << "********************************************************************************" << std::endl;

	for (std::size_t i = 0; i < _registered_benchmarks.size(); i++)
	{
		std::cout << _registered_benchmarks[i]->getHelp();
		std::cout << std::endl;
	}

	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << "Benchmarks (END)" << std::endl;
	o_ostream << "********************************************************************************" << std::endl;

}


BenchmarksCombined::~BenchmarksCombined()
{
	clear();
}

}}
