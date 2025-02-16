/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONSPHERE2D_TIME_PDEADVECTIONSPHERE2DTS_BASEINTERFACE_HPP
#define PROGRAMS_PDE_ADVECTIONSPHERE2D_TIME_PDEADVECTIONSPHERE2DTS_BASEINTERFACE_HPP

#include <programs/PDE_AdvectionSphere2D/time/ShackPDEAdvectionSphere2DTimeDiscretization.hpp>
#include <sweet/Data/Sphere2D/Shack.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/SemiLagrangian/Shack.hpp>
#include <limits>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>

#include "../benchmarks/ShackPDEAdvectionSphere2DBenchmarks.hpp"


class PDEAdvectionSphere2DTS_BaseInterface
{
public:
	sweet::Error::Base error;

	/*
	 * These are just some default shacks we provide to each time stepping method
	 */
	sweet::Shacks::Dictionary *shackDict;
	sweet::TimeTree::Shack *shackTimestepControl;
	sweet::Data::Sphere2D::Shack *shackSphere2DDataOps;
	ShackPDEAdvectionSphere2DTimeDiscretization *shackPDEAdvectionTimeDisc;
	ShackPDEAdvectionSphere2DBenchmarks *shackPDEAdvBenchmark;
	sweet::SemiLagrangian::Shack *shackSemiLagrangian;

	sweet::Data::Sphere2D::Operators *ops;

	PDEAdvectionSphere2DTS_BaseInterface()	:
		shackDict(nullptr),
		shackTimestepControl(nullptr),
		shackSphere2DDataOps(nullptr),
		shackPDEAdvectionTimeDisc(nullptr),
		shackPDEAdvBenchmark(nullptr),
		ops(nullptr)
	{

	}

	virtual bool shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackTimestepControl = io_shackDict->getAutoRegistration<sweet::TimeTree::Shack>();
		shackSphere2DDataOps = io_shackDict->getAutoRegistration<sweet::Data::Sphere2D::Shack>();
		shackPDEAdvectionTimeDisc = io_shackDict->getAutoRegistration<ShackPDEAdvectionSphere2DTimeDiscretization>();
		shackPDEAdvBenchmark = io_shackDict->getAutoRegistration<ShackPDEAdvectionSphere2DBenchmarks>();
		shackSemiLagrangian = io_shackDict->getAutoRegistration<sweet::SemiLagrangian::Shack>();

		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}

	virtual bool setup(
			sweet::Data::Sphere2D::Operators *io_ops
	)
	{
		ops = io_ops;
		return true;
	}

	/*
	 * Timestepping interface used by main timestepping loop
	 */
	virtual void runTimestep(
			std::vector<sweet::Data::Sphere2D::DataSpectral> &io_prognostic_fields,	//!< prognostic variables
			sweet::Data::Sphere2D::DataGrid &io_u,
			sweet::Data::Sphere2D::DataGrid &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp
	) = 0;

	virtual bool testImplementsTimesteppingMethod(
			const std::string &i_timestepping_method
		) = 0;

	virtual std::string getStringId() = 0;

	virtual void printImplementedTimesteppingMethods(
			std::ostream &o_ostream = std::cout,
			const std::string &i_prefix = ""
		) = 0;

	virtual ~PDEAdvectionSphere2DTS_BaseInterface()
	{
	}
};

#endif
