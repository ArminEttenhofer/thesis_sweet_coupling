/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONCART2D_PDEADVECTIONCART2DBENCHMARKSCOMBINED_HPP
#define PROGRAMS_PDE_ADVECTIONCART2D_PDEADVECTIONCART2DBENCHMARKSCOMBINED_HPP

#include <programs/PDE_AdvectionCart2D/benchmarks/PDEAdvectionCart2DBench_Cylinder.hpp>
#include <programs/PDE_AdvectionCart2D/benchmarks/PDEAdvectionCart2DBench_GaussianBump.hpp>
#include <programs/PDE_AdvectionCart2D/benchmarks/PDEAdvectionCart2DBench_GaussianBumpAdvection.hpp>
#include <programs/PDE_AdvectionCart2D/benchmarks/PDEAdvectionCart2DBench_RadialGaussianBump.hpp>
#include <programs/PDE_AdvectionCart2D/benchmarks/ShackPDEAdvectionCart2DBenchmarks.hpp>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <iostream>
#include <sweet/Error/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>




class PDEAdvectionCart2DBenchmarksCombined
{
public:
	sweet::Error::Base error;
	ShackPDEAdvectionCart2DBenchmarks *shackBenchmarks;
	sweet::Data::Cart2D::Operators *op;

	PDEAdvectionCart2DBenchmarksCombined()	:
		shackBenchmarks(nullptr),
		op(nullptr)
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
		shackBenchmarks = io_shackDict.getAutoRegistration<ShackPDEAdvectionCart2DBenchmarks>();

		return error.forwardWithPositiveReturn(io_shackDict.error);
	}

	void clear()
	{
		shackBenchmarks = nullptr;
	}


public:

public:
	bool setupInitialConditions(
			sweet::Data::Cart2D::DataSpectral &o_h_pert,
			sweet::Data::Cart2D::DataSpectral &o_u,
			sweet::Data::Cart2D::DataSpectral &o_v,
			sweet::Data::Cart2D::Operators &io_op,				//!< Make this IO, since changes in the simulation parameters might require to also update the operators
			sweet::Shacks::Dictionary &io_shackDict
	)
	{
		if (!shackBenchmarks->validateNonzeroAdvection())
			return error.forwardWithPositiveReturn(shackBenchmarks->error);

		if (shackBenchmarks->benchmark_name == "")
			return error.set("SWECart2DBenchmarksCombined: Benchmark name not given, use --benchmark-name=[name]");


		if (shackBenchmarks->benchmark_name == "gaussian_bump" || shackBenchmarks->benchmark_name == "gaussian_bump_phi_pint")
		{
			PDEAdvectionCart2DBenchGaussianBump gaussian_bump;

			gaussian_bump.shackRegistration(&io_shackDict);

			gaussian_bump.setup(op);

			gaussian_bump.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (shackBenchmarks->benchmark_name == "gaussian_bump_advection")
		{
			PDEAdvectionCart2DBenchGaussianBumpAdvection gaussian_bump_adv;

			gaussian_bump_adv.shackRegistration(&io_shackDict);

			gaussian_bump_adv.setup(op);

			gaussian_bump_adv.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (
				shackBenchmarks->benchmark_name == "benchmark_id_0" ||
				shackBenchmarks->benchmark_name == "cylinder"
		)
		{

			PDEAdvectionCart2DBenchCylinder gaussian_bump_cylinder;

			gaussian_bump_cylinder.shackRegistration(&io_shackDict);

			gaussian_bump_cylinder.setup(op);

			gaussian_bump_cylinder.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (
				shackBenchmarks->benchmark_name == "benchmark_id_1" ||
				shackBenchmarks->benchmark_name == "radial_gaussian_bump"
		)
		{

			PDEAdvectionCart2DBenchRadialGaussianBump gaussian_radial_gaussian_bump;

			gaussian_radial_gaussian_bump.shackRegistration(&io_shackDict);

			gaussian_radial_gaussian_bump.setup(op);

			gaussian_radial_gaussian_bump.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}

		shackBenchmarks->printProgramArguments();

		return error.set(std::string("Benchmark '")+shackBenchmarks->benchmark_name+ "' not found (or not available)");
	}
};



#endif
