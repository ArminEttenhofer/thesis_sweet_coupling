/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_BENCHMARKSCOMBINED_HPP
#define PROGRAMS_PDE_SWESPHERE2D_BENCHMARKSCOMBINED_HPP

#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include "Benchmarks/BaseInterface.hpp"
#include "Benchmarks/Shack.hpp"

#if !SWEET_USE_SPHERE2D_SPECTRAL_SPACE && SWEET_PARAREAL
	#error "SWEET_USE_SPHERE2D_SPECTRAL_SPACE not enabled"
#endif



#include <iostream>
#include <sweet/Error/Base.hpp>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Shacks/Dictionary.hpp>

#include "Benchmarks/Shack.hpp"
#include "Benchmarks/BaseInterface.hpp"
#include <sweet/TimeTree/Shack.hpp>

#include <sweet/Data/Sphere2D/Sphere2D.hpp>


namespace PDE_SWESphere2D {
namespace Benchmarks {

class BenchmarksCombined
{
	std::vector<BaseInterface*> _registered_benchmarks;

	sweet::Shacks::Dictionary *shackDict;
	Shack *shackBenchmarks;

public:
	sweet::Error::Base error;

	// Sphere2D operators
	sweet::Data::Sphere2D::Operators *ops;

	BaseInterface *benchmark = nullptr;

	BenchmarksCombined();

	bool setup_1_registerAllBenchmark();

	bool setup_2_shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
	);

	bool setup_3_benchmarkDetection(
			const std::string &i_benchmark_name = ""
	);

	bool setup_4_benchmarkSetup_1_withoutOps();

	bool setup_5_benchmarkSetup_2_withOps(
			sweet::Data::Sphere2D::Operators *io_ops
	);

	void clear_3_benchmarkDetection();
	void clear();

private:
	void _benchmarksFreeAll(
			BaseInterface *skip_this = nullptr
	);

public:
	void printAvailableBenchmarks(
		std::ostream &o_ostream = std::cout,
		const std::string &i_prefix = ""
	);

public:
	~BenchmarksCombined();
};

}}

#endif
