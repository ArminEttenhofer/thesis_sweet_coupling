#ifndef PROGRAMS_ODE_GENERIC_DE_BASE_HPP
#define PROGRAMS_ODE_GENERIC_DE_BASE_HPP

#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Data/GenericContainer/ConfigBase.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>

namespace ODE_Generic {

/*!
 * Base class to abstract away the different ODEs
 */
class DE_Base
{
public:
	sweet::Error::Base error;

	DE_Base()
	{
	}

	virtual ~DE_Base()
	{
	}

	virtual
	bool setup(
			sweet::Shacks::ProgramArgumentsDictionary *i_progArgShackDict,
			const std::string &i_timesteppingMethod
	)
	{
		return error.set("TODO: setup(...) needs to be overridden");
	}

	virtual
	bool setTimeStepSize(
			double i_dt
	)
	{
		return error.set("TODO: setTimeStepSize(...) needs to be overridden");
	}

	virtual
	bool getInitialState(
			sweet::Data::GenericContainer::Base &o_U
	)
	{
		return error.set("TODO: getInitialState(...) needs to be overridden");
	}

	virtual
	bool runTimestep(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		return error.set("TODO: runTimestep(...) needs to be overridden");
	}

	virtual
	bool runFileOutput(
			const sweet::Data::GenericContainer::Base &i_data
	)
	{
		return error.set("TODO: runFileOutput(...) needs to be overridden");
	}

	virtual
	bool clear()
	{
		return error.set("TODO: clear() needs to be overridden");
	}


	virtual
	bool print(
			const sweet::Data::GenericContainer::Base &i_data,
			std::string &i_type,
			std::iostream &io_iostream
	)
	{
		return error.set("TODO: print(...) needs to be overridden");
		return true;
	}

	virtual
	double computeError(
			const sweet::Data::GenericContainer::Base &i_data,
			double i_timeStamp,
			const std::string &i_errorNorm
	)
	{
		return error.set("TODO: computeError(...) needs to be overridden");
	}

	virtual
	sweet::Data::GenericContainer::Base* getNewDataContainerInstance(int i_id = sweet::Data::GenericContainer::Base::DATA_SIMULATION)
	{
		error.set("TODO: getNewDataContainerInstance(...) needs to be overridden");
		return nullptr;
	}
};

}

#endif
