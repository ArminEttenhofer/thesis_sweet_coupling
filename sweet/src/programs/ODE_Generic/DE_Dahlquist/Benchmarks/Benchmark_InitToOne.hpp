/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_BENCHMARKS_BENCHMARK_INITTOONE_HPP
#define PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_BENCHMARKS_BENCHMARK_INITTOONE_HPP

#include <sweet/Data/GenericContainer/CastHelper.hpp>
#include <sweet/Shacks/Dictionary.hpp>

#include "../DataContainer/Config.hpp"
#include "Base.hpp"
#include "../Shack.hpp"
#include "../DataContainer/Simulation.hpp"

namespace ODE_Generic {
namespace DE_Dahlquist {
namespace Benchmarks {


class Benchmark_InitToOne	:
		public Base,
		public sweet::Data::GenericContainer::CastHelper<ODE_Generic::DE_Dahlquist::DataContainer::Simulation, ODE_Generic::DE_Dahlquist::DataContainer::Config>
{
public:
	DE_Dahlquist::Shack *shackDahlquist;

	Benchmark_InitToOne()	:
		shackDahlquist(nullptr)
	{
	}

	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override
	{
		shackDahlquist = io_shackDict->getAutoRegistration<DE_Dahlquist::Shack>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

		Base::shackRegistration(io_shackDict);
		return true;
	}

	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		) override
	{
		benchmark_name = i_benchmark_name;

		return
			i_benchmark_name == "initone"	||
			false
		;
	}


	void setup_1_shackData() override
	{
	}

	void clear() override
	{
	}

	std::string printHelp() override
	{
		std::ostringstream stream;
		stream << "  'initone': Set U=1" << std::endl;
		return stream.str();
	}


	/*!
	 * Simply set U to 1.0
	 */
	void getInitialState(
		sweet::Data::GenericContainer::Base &o_U_
	) override
	{
		cast(o_U_).U = 1.0;
	}

	virtual void getReferenceState(
		sweet::Data::GenericContainer::Base &o_U_,
		double i_timeStamp
	) override
	{
		SWEET_ASSERT(shackDahlquist != nullptr);

		getInitialState(o_U_);

		DataContainer::Simulation &o_U = cast(o_U_);
		std::complex<double> sumLambda = shackDahlquist->lambda1 + shackDahlquist->lambda2 + shackDahlquist->lambda3;

		o_U.U = std::exp(i_timeStamp*sumLambda)*o_U.U;
	}
};

}}}

#endif
