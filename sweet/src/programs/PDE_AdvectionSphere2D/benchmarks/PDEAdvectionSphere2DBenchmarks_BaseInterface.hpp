/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONSPHERE2D_BENCHMARKS_PDEADVECTIONSPHERE2DBENCHMARKS_BASEINTERFACE_HPP
#define PROGRAMS_PDE_ADVECTIONSPHERE2D_BENCHMARKS_PDEADVECTIONSPHERE2DBENCHMARKS_BASEINTERFACE_HPP

#include <programs/PDE_AdvectionSphere2D/time/ShackPDEAdvectionSphere2DTimeDiscretization.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Data/Vector/Vector.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>

#include "../benchmarks/ShackPDEAdvectionSphere2DBenchmarks.hpp"


class PDEAdvectionSphere2DBenchmarks_BaseInterface
{
public:
	sweet::Error::Base error;

	/*
	 * These are just some default shacks we provide to each time stepping method
	 */
	sweet::Shacks::Dictionary *shackDict;
	sweet::TimeTree::Shack *shackTimestepControl;
	sweet::Data::Sphere2D::Shack *shackSphere2DDataOps;
	ShackPDEAdvectionSphere2DTimeDiscretization *shackPDEAdvTimeDisc;
	ShackPDEAdvectionSphere2DBenchmarks *shackPDEAdvBenchmark;

	sweet::Data::Sphere2D::Operators *ops;

	PDEAdvectionSphere2DBenchmarks_BaseInterface()	:
		shackDict(nullptr),
		shackTimestepControl(nullptr),
		shackSphere2DDataOps(nullptr),
		shackPDEAdvTimeDisc(nullptr),
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
		shackPDEAdvTimeDisc = io_shackDict->getAutoRegistration<ShackPDEAdvectionSphere2DTimeDiscretization>();
		shackPDEAdvBenchmark = io_shackDict->getAutoRegistration<ShackPDEAdvectionSphere2DBenchmarks>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}

public:
	virtual void setup_1_shackData(
	) = 0;

public:
	virtual void setup_2_withOps(
			sweet::Data::Sphere2D::Operators *io_ops
	) = 0;

	virtual bool implements_benchmark(
			const std::string &i_benchmark_name
		) = 0;


	virtual std::string printHelp() = 0;

	virtual void getInitialState(
		std::vector<sweet::Data::Sphere2D::DataSpectral> &o_phi,
		sweet::Data::Sphere2D::DataGrid &o_u,
		sweet::Data::Sphere2D::DataGrid &o_v
	) = 0;


#if 0
	virtual void get_reference_state(
		std::vector<sweet::Data::Sphere2D::DataSpectral*> &o_phi_pert,
		double i_timestamp
	)
	{
		SWEETErrorFatal("'get_reference_state' for multiple prognostic variables not implemented for this benchmark");
	}


	virtual void get_varying_velocities(
		sweet::Data::Sphere2D::DataGrid &o_u,
		sweet::Data::Sphere2D::DataGrid &o_v,
		double i_timestamp
	)
	{
		SWEETErrorFatal("'get_varying_velocities' not implemented for this benchmark");
	}


	virtual bool has_time_varying_state()
	{
		return false;
	}
#endif


	/*
	 * Return number of prognostic fields to be used
	 */
	virtual int getNumPrognosticFields()
	{
		return 1;
	}


	virtual ~PDEAdvectionSphere2DBenchmarks_BaseInterface()
	{
	}
};

#endif
