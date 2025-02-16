/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_BASEINTERFACE_HPP
#define PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_BASEINTERFACE_HPP

#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Data/Vector/Vector.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/Parallelization/Shack.hpp>

#include "Shack.hpp"
#include "../Shack.hpp"
#include <sweet/TimeTree/Shack.hpp>

#include "../TimeOld/ShackTimeDiscretization.hpp"

namespace PDE_SWESphere2D {
namespace Benchmarks {

class BaseInterface
{
public:
	sweet::Error::Base error;

	/*
	 * These are just some default shacks we provide to each time stepping method
	 */
	sweet::Shacks::Dictionary *shackDict;
	sweet::TimeTree::Shack *shackTimestepControl;
	sweet::Data::Sphere2D::Shack *shackSphere2DDataOps;
	sweet::Parallelization::Shack *shackParallelization;
	ShackTimeDiscretization *shackPDESWETimeDisc;
	Benchmarks::Shack *shackPDESWEBenchmarks;
	PDE_SWESphere2D::Shack *shackPDESWESphere2D;

	sweet::Data::Sphere2D::Operators *ops;

	BaseInterface()	:
		shackDict(nullptr),
		shackTimestepControl(nullptr),
		shackSphere2DDataOps(nullptr),
		shackPDESWETimeDisc(nullptr),
		shackPDESWEBenchmarks(nullptr),
		shackPDESWESphere2D(nullptr),
		ops(nullptr)
	{
	}

	virtual bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackDict->getAutoRegistration(&shackTimestepControl);
		shackDict->getAutoRegistration(&shackSphere2DDataOps);
		shackDict->getAutoRegistration(&shackPDESWETimeDisc);
		shackDict->getAutoRegistration(&shackPDESWEBenchmarks);
		shackDict->getAutoRegistration(&shackPDESWESphere2D);
		shackDict->getAutoRegistration(&shackParallelization);

		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*shackDict);

		return true;
	}

public:
	virtual void setup_1_shackData() = 0;

public:
	virtual void setup_2_withOps(
			sweet::Data::Sphere2D::Operators *io_ops
	) = 0;

	virtual bool implements_benchmark(
			const std::string &i_benchmark_name
		) = 0;


	virtual std::string getHelp() = 0;

	virtual void getInitialState(
		sweet::Data::Sphere2D::DataSpectral &o_phi,
		sweet::Data::Sphere2D::DataSpectral &o_vrt,
		sweet::Data::Sphere2D::DataSpectral &o_div
	) = 0;

	virtual void setup_topography(
	) = 0;

	virtual void getReferenceState(
		sweet::Data::Sphere2D::DataSpectral &o_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &o_vrt,
		sweet::Data::Sphere2D::DataSpectral &o_div,
		double i_timestamp
	)
	{
		SWEETErrorFatal("Not implemented for this benchmark");
	}

	virtual bool has_time_varying_state()
	{
		return false;
	}

	virtual void clear() = 0;

	virtual ~BaseInterface()
	{
	}
};

}}

#endif
