/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWECART2D_BENCHMARKS_PDESWECART2DBENCH_BASEINTERFACE_HPP
#define PROGRAMS_PDE_SWECART2D_BENCHMARKS_PDESWECART2DBENCH_BASEINTERFACE_HPP


#include <programs/PDE_SWECart2D/Benchmarks/Shack_Polvani.hpp>
#include <programs/PDE_SWECart2D/Benchmarks/Shack.hpp>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>

#include "../Shack.hpp"

namespace PDE_SWECart2D {
namespace Benchmarks {

class BaseInterface
{
public:
	sweet::Error::Base error;

	sweet::Shacks::Dictionary *shackDict;
	sweet::Data::Cart2D::Shack *shackCart2DDataOps;
	sweet::TimeTree::Shack *shackTimestepControl;
	PDE_SWECart2D::Benchmarks::Shack *shackBenchmarks;
	PDE_SWECart2D::Shack *shackPDESWECart2D;
	PDE_SWECart2D::Benchmarks::Shack_Polvani *shackPDESWECart2DBench_PolvaniBench;

	sweet::Data::Cart2D::Operators *ops;
	sweet::Data::Cart2D::Config *cart2DDataConfig;


	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		////shackDict->getAutoRegistration(&shackCart2DDataOps);
		////shackDict->getAutoRegistration(&shackTimestepControl);
		////shackDict->getAutoRegistration(&shackBenchmarks);
		////shackDict->getAutoRegistration(&shackPDESWECart2D);
		////shackDict->getAutoRegistration(&shackPDESWECart2DBench_PolvaniBench);

		shackCart2DDataOps = io_shackDict->getAutoRegistration<sweet::Data::Cart2D::Shack>();
		shackTimestepControl = io_shackDict->getAutoRegistration<sweet::TimeTree::Shack>();
		shackBenchmarks = io_shackDict->getAutoRegistration<PDE_SWECart2D::Benchmarks::Shack>();
		shackPDESWECart2D = io_shackDict->getAutoRegistration<PDE_SWECart2D::Shack>();
		shackPDESWECart2DBench_PolvaniBench = io_shackDict->getAutoRegistration<PDE_SWECart2D::Benchmarks::Shack_Polvani>();

		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(*io_shackDict);
	}

	virtual bool setup(
		sweet::Data::Cart2D::Operators *io_ops,
		sweet::Data::Cart2D::Config *i_cart2DDataConfig
	)
	{
		ops = io_ops;
		cart2DDataConfig = i_cart2DDataConfig;

		return true;
	}

	virtual bool setupBenchmark(
			sweet::Data::Cart2D::DataSpectral &o_h_pert,
			sweet::Data::Cart2D::DataSpectral &o_u,
			sweet::Data::Cart2D::DataSpectral &o_v
	) = 0;
};

}}

#endif
