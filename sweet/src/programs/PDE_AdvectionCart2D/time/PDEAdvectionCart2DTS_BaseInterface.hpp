/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONCART2D_TIME_PDEADVECTIONCART2DTS_BASEINTERFACE_HPP
#define PROGRAMS_PDE_ADVECTIONCART2D_TIME_PDEADVECTIONCART2DTS_BASEINTERFACE_HPP

#include <programs/PDE_AdvectionCart2D/benchmarks/ShackPDEAdvectionCart2DBenchmarks.hpp>
#include <programs/PDE_AdvectionCart2D/time/ShackPDEAdvectionCart2DTimeDiscretization.hpp>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <limits>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>



class PDEAdvectionCart2DTS_BaseInterface
{
public:
	sweet::Error::Base error;

	/*
	 * These are just some default shacks we provide to each time stepping method
	 */
	sweet::Shacks::Dictionary *shackDict;
	sweet::TimeTree::Shack *shackTimestepControl;
	sweet::Data::Cart2D::Shack *shackCart2DDataOps;
	ShackPDEAdvectionCart2DTimeDiscretization *shackPDEAdvTimeDisc;
	ShackPDEAdvectionCart2DBenchmarks *shackPDEAdvBenchmark;

	sweet::Data::Cart2D::Operators *ops;

	PDEAdvectionCart2DTS_BaseInterface()	:
		shackDict(nullptr),
		shackTimestepControl(nullptr),
		shackCart2DDataOps(nullptr),
		shackPDEAdvTimeDisc(nullptr),
		shackPDEAdvBenchmark(nullptr),
		ops(nullptr)
	{

	}

	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackTimestepControl = io_shackDict->getAutoRegistration<sweet::TimeTree::Shack>();
		shackCart2DDataOps = io_shackDict->getAutoRegistration<sweet::Data::Cart2D::Shack>();
		shackPDEAdvTimeDisc = io_shackDict->getAutoRegistration<ShackPDEAdvectionCart2DTimeDiscretization>();
		shackPDEAdvBenchmark = io_shackDict->getAutoRegistration<ShackPDEAdvectionCart2DBenchmarks>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}

	virtual bool setup(
			sweet::Data::Cart2D::Operators *io_ops
	)
	{
		ops = io_ops;
		return true;
	}


public:
	virtual void runTimestep(
			sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	) = 0;

	~PDEAdvectionCart2DTS_BaseInterface() {}
};



#endif
