/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_BENCHMARKS_BASE_HPP
#define PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_BENCHMARKS_BASE_HPP

#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>

#include "Shack.hpp"

namespace ODE_Generic {
namespace DE_Dahlquist {
namespace Benchmarks {

class Base
{
public:
	sweet::Error::Base error;

	/*
	 * These are just some default shacks we provide to each time stepping method
	 */
	sweet::Shacks::Dictionary *shackDict;
	sweet::TimeTree::Shack *shackTimestepControl;
	Shack *shackODEScalarBenchmarks;
	Shack *shackODEScalar;

	Base()	:
		shackDict(nullptr),
		shackTimestepControl(nullptr),
		shackODEScalarBenchmarks(nullptr),
		shackODEScalar(nullptr)
	{
	}

	virtual bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackTimestepControl = shackDict->getAutoRegistration<sweet::TimeTree::Shack>();
		shackODEScalarBenchmarks = shackDict->getAutoRegistration<Shack>();
		shackODEScalar = shackDict->getAutoRegistration<Shack>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*shackDict);

		return true;
	}

public:
	virtual void setup_1_shackData() = 0;

public:
	virtual bool implements_benchmark(
			const std::string &i_benchmark_name
		) = 0;


	virtual std::string printHelp() = 0;

	virtual void getInitialState(
		sweet::Data::GenericContainer::Base &o_U
	) = 0;


	virtual void getReferenceState(
		sweet::Data::GenericContainer::Base &o_U,
		double i_timeStamp
	)
	{
		SWEETErrorFatal("Not implemented for this benchmark");
	}

	virtual bool has_time_varying_state()
	{
		return false;
	}

	virtual void clear() = 0;

	virtual ~Base()
	{
	}
};

}}}

#endif
