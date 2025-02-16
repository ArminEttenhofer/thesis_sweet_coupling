/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <programs/PDE_AdvectionSphere2D/benchmarks/PDEAdvectionSphere2DBenchmark_advection_vector_3d_normal_vectors.hpp>
#include <programs/PDE_AdvectionSphere2D/benchmarks/PDEAdvectionSphere2DBenchmark_advection_vector_uv_gauss_bumps.hpp>
#include <programs/PDE_AdvectionSphere2D/benchmarks/PDEAdvectionSphere2DBenchmark_advection_vector_uv_velocities.hpp>
#include <programs/PDE_AdvectionSphere2D/benchmarks/PDEAdvectionSphere2DBenchmark_nair_lauritzen_sl.hpp>
#include <programs/PDE_AdvectionSphere2D/benchmarks/PDEAdvectionSphere2DBenchmark_williamson_1_advection_cos_bell.hpp>
#include <programs/PDE_AdvectionSphere2D/benchmarks/PDEAdvectionSphere2DBenchmark_williamson_1_advection_gauss_bump.hpp>
#include <programs/PDE_AdvectionSphere2D/benchmarks/PDEAdvectionSphere2DBenchmark_zero.hpp>
#include <programs/PDE_AdvectionSphere2D/PDEAdvectionSphere2DBenchmarksCombined.hpp>
#include <programs/PDE_AdvectionSphere2D/time/PDEAdvectionSphere2DTS_BaseInterface.hpp>



PDEAdvectionSphere2DBenchmarksCombined::PDEAdvectionSphere2DBenchmarksCombined()	:
	shackDict(nullptr),
	ops(nullptr),
	benchmark(nullptr)
{
}


bool PDEAdvectionSphere2DBenchmarksCombined::setup_1_registerAllBenchmark()
{
	_registered_benchmarks.push_back(static_cast<PDEAdvectionSphere2DBenchmarks_BaseInterface*>(new PDEAdvectionSphere2DBenchmark_zero));
	_registered_benchmarks.push_back(static_cast<PDEAdvectionSphere2DBenchmarks_BaseInterface*>(new PDEAdvectionSphere2DBenchmark_nair_lauritzen_sl));
	_registered_benchmarks.push_back(static_cast<PDEAdvectionSphere2DBenchmarks_BaseInterface*>(new PDEAdvectionSphere2DBenchmark_advection_vector_uv_velocities));
	_registered_benchmarks.push_back(static_cast<PDEAdvectionSphere2DBenchmarks_BaseInterface*>(new PDEAdvectionSphere2DBenchmark_advection_vector_uv_gauss_bumps));
	_registered_benchmarks.push_back(static_cast<PDEAdvectionSphere2DBenchmarks_BaseInterface*>(new PDEAdvectionSphere2DBenchmark_advection_vector_3d_normal_vectors));
	_registered_benchmarks.push_back(static_cast<PDEAdvectionSphere2DBenchmarks_BaseInterface*>(new PDEAdvectionSphere2DBenchmark_williamson_1_advection_cos_bell));
	_registered_benchmarks.push_back(static_cast<PDEAdvectionSphere2DBenchmarks_BaseInterface*>(new PDEAdvectionSphere2DBenchmark_williamson_1_advection_gauss_bump));

	return true;
}



bool PDEAdvectionSphere2DBenchmarksCombined::setup_2_shackRegistration(
	sweet::Shacks::Dictionary *io_shackDict
)
{
	shackDict = io_shackDict;

	shackBenchmarks = shackDict->getAutoRegistration<ShackPDEAdvectionSphere2DBenchmarks>();

	for (std::size_t i = 0; i < _registered_benchmarks.size(); i++)
	{
		_registered_benchmarks[i]->shackRegistration(io_shackDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_registered_benchmarks[i]);
	}

	return true;
}



bool PDEAdvectionSphere2DBenchmarksCombined::setup_3_benchmarkDetection()
{
	SWEET_ASSERT(benchmark == nullptr);

	if (shackBenchmarks->benchmark_name == "")
	{
		printAvailableBenchmarks();
		return error.set("Please choose benchmark with --benchmark-name=...");
	}

	/*
	 * Find right one
	 */
	for (std::size_t i = 0; i < _registered_benchmarks.size(); i++)
	{
		PDEAdvectionSphere2DBenchmarks_BaseInterface *ts = _registered_benchmarks[i];

		if (ts->implements_benchmark(shackBenchmarks->benchmark_name))
		{
			if (benchmark != nullptr)
			{
				//std::cout << "Processing " << i+1 << "th element" << std::endl;
				error.set("Duplicate implementation for benchmark "+shackBenchmarks->benchmark_name);
			}

			//std::cout << "Benchmark detection: found match with benchmark id " << i+1 << std::endl;
			benchmark = ts;
		}
	}

	if (benchmark == nullptr)
	{
		printAvailableBenchmarks();
		return error.set("No valid --benchmark-name '"+shackBenchmarks->benchmark_name+"' provided");
	}

	return true;
}


bool PDEAdvectionSphere2DBenchmarksCombined::setup_4_benchmarkSetup_1_withoutOps()
{
	SWEET_ASSERT(benchmark != nullptr);

	benchmark->setup_1_shackData();
	ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(*benchmark);

	return true;
}


bool PDEAdvectionSphere2DBenchmarksCombined::setup_5_benchmarkSetup_2_withOps(
		sweet::Data::Sphere2D::Operators* io_ops
)
{
	ops = io_ops;

	benchmark->setup_2_withOps(ops);
	ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(*benchmark);
}

void PDEAdvectionSphere2DBenchmarksCombined::clear()
{
	benchmark = nullptr;
	shackDict = nullptr;
	ops = nullptr;

	_benchmarksFreeAll();
}



void PDEAdvectionSphere2DBenchmarksCombined::_benchmarksFreeAll(
		PDEAdvectionSphere2DBenchmarks_BaseInterface *skip_this
)
{
	for (std::size_t i = 0; i < _registered_benchmarks.size(); i++)
	{
		PDEAdvectionSphere2DBenchmarks_BaseInterface *ts = _registered_benchmarks[i];

		if (ts == skip_this)
			continue;

		delete ts;
	}

	_registered_benchmarks.clear();
}


void PDEAdvectionSphere2DBenchmarksCombined::printAvailableBenchmarks(
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


PDEAdvectionSphere2DBenchmarksCombined::~PDEAdvectionSphere2DBenchmarksCombined()
{
	_benchmarksFreeAll();
}
