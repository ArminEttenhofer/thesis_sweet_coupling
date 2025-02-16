/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONCART2D_BENCHMARKS_PDEADVECTIONCART2DBENCH_BASEINTERFACE_HPP
#define PROGRAMS_PDE_ADVECTIONCART2D_BENCHMARKS_PDEADVECTIONCART2DBENCH_BASEINTERFACE_HPP


#include <programs/PDE_AdvectionCart2D/benchmarks/ShackPDEAdvectionCart2DBenchmarks.hpp>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>



class PDEAdvectionCart2DBench_BaseInterface
{
public:
	sweet::Error::Base error;

	sweet::Shacks::Dictionary *shackDict;
	sweet::Data::Cart2D::Shack *shackCart2DDataOps;
	sweet::TimeTree::Shack *shackTimestepControl;
	ShackPDEAdvectionCart2DBenchmarks *shackBenchmarks;

	sweet::Data::Cart2D::Operators *ops;


	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackCart2DDataOps = io_shackDict->getAutoRegistration<sweet::Data::Cart2D::Shack>();
		shackTimestepControl = io_shackDict->getAutoRegistration<sweet::TimeTree::Shack>();
		shackBenchmarks = io_shackDict->getAutoRegistration<ShackPDEAdvectionCart2DBenchmarks>();

		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(*io_shackDict);
	}

	virtual bool setup(
		sweet::Data::Cart2D::Operators *io_ops
	)
	{
		ops = io_ops;

		return true;
	}

	virtual bool setupBenchmark(
			sweet::Data::Cart2D::DataSpectral &o_h_pert,
			sweet::Data::Cart2D::DataSpectral &o_u,
			sweet::Data::Cart2D::DataSpectral &o_v
	) = 0;
};


#endif
