/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONSPHERE2D_PDEADVECTIONSPHERE2DBENCHMARKSCOMBINED_HPP
#define PROGRAMS_PDE_ADVECTIONSPHERE2D_PDEADVECTIONSPHERE2DBENCHMARKSCOMBINED_HPP


#include <programs/PDE_AdvectionSphere2D/benchmarks/PDEAdvectionSphere2DBenchmarks_BaseInterface.hpp>
#include <programs/PDE_AdvectionSphere2D/benchmarks/ShackPDEAdvectionSphere2DBenchmarks.hpp>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <iostream>
#include <sweet/Error/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>



class PDEAdvectionSphere2DBenchmarksCombined
{
	std::vector<PDEAdvectionSphere2DBenchmarks_BaseInterface*> _registered_benchmarks;

	sweet::Shacks::Dictionary *shackDict;
	ShackPDEAdvectionSphere2DBenchmarks *shackBenchmarks;

public:
	sweet::Error::Base error;

	// Sphere2D operators
	sweet::Data::Sphere2D::Operators *ops;

	PDEAdvectionSphere2DBenchmarks_BaseInterface *benchmark = nullptr;

	PDEAdvectionSphere2DBenchmarksCombined();

	bool setup_1_registerAllBenchmark();

	bool setup_2_shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
	);

	bool setup_3_benchmarkDetection();

	bool setup_4_benchmarkSetup_1_withoutOps();

	bool setup_5_benchmarkSetup_2_withOps(
			sweet::Data::Sphere2D::Operators *io_ops
	);

	void clear();

private:
	void _benchmarksFreeAll(
			PDEAdvectionSphere2DBenchmarks_BaseInterface *skip_this = nullptr
	);

public:
	void printAvailableBenchmarks(
		std::ostream &o_ostream = std::cout,
		const std::string &i_prefix = ""
	);

public:
	~PDEAdvectionSphere2DBenchmarksCombined();
};


#endif
