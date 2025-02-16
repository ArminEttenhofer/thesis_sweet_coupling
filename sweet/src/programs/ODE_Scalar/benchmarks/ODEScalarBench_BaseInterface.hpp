/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef PROGRAMS_ODE_SCALAR_BENCHMARKS_ODESCALARBENCH_BASEINTERFACE_HPP
#define PROGRAMS_ODE_SCALAR_BENCHMARKS_ODESCALARBENCH_BASEINTERFACE_HPP


#include <programs/ODE_Scalar/benchmarks/ShackODEScalarBenchmarks.hpp>
#include <programs/ODE_Scalar/ShackODEScalar.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>


class ODEScalarBench_BaseInterface
{
public:
	sweet::Error::Base error;

	sweet::Shacks::Dictionary *shackDict;
	sweet::TimeTree::Shack *shackTimestepControl;
	ShackODEScalar *shackODEScalar;
	ShackODEScalarBenchmarks *shackBenchmarks;

	ODEScalarBench_BaseInterface() :
		shackDict(nullptr),
		shackTimestepControl(nullptr),
		shackODEScalar(nullptr),
		shackBenchmarks(nullptr)
	{
	}

	virtual bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackTimestepControl = io_shackDict->getAutoRegistration<sweet::TimeTree::Shack>();
		shackODEScalar = io_shackDict->getAutoRegistration<ShackODEScalar>();
		shackBenchmarks = io_shackDict->getAutoRegistration<ShackODEScalarBenchmarks>();

		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}

	virtual bool setup(
	)
	{
		return true;
	}

	virtual bool setupBenchmark(
			double &o_u
	) = 0;
};


#endif

