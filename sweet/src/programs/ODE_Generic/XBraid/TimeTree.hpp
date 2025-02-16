#ifndef PROGRAMS_ODE_GENERIC_XBRAID_TIMETREE_HPP
#define PROGRAMS_ODE_GENERIC_XBRAID_TIMETREE_HPP

#include <sweet/Error/Base.hpp>
#include <sweet/XBraid/TimeTree.hpp>
#include <sweet/XBraid/Vector.hpp>
#include <programs/ODE_Generic/DE_Dahlquist/TimeTree/TimeTreeIR.hpp>
#include <programs/ODE_Generic/DE_Dahlquist/DataContainer/Simulation.hpp>
#include <programs/ODE_Generic/XBraid/DataContainer.hpp>


namespace ODE_Generic {
namespace XBraid {

class TimeTree :
	public sweet::XBraid::TimeTree
{

public:

	sweet::Error::Base error;

	ODE_Generic::TimeTree::Base* tsm = nullptr;
	sweet::Data::GenericContainer::ConfigBase* config = nullptr;

public:

	TimeTree()
	{
		/////tsm = new ODE_Generic::TimeTree::TimeTree;
	}

public:
	~TimeTree()
	{
		clear();
	}

public:
	void clear()
	{
		if (tsm)
		{
			tsm->clear();
			tsm = nullptr;
		}

		if (config)
		{
			delete config;
			config = nullptr;
		}
	}

public:
	void setup(
			std::string &i_timestepping_method,
			sweet::Shacks::Dictionary* i_shackDict,
			ODE_Generic::XBraid::DataContainer* i_u,
			double i_dt
	)
	{
		tsm->setup_1_registerAllTimesteppers();
		tsm->setup_2_shackRegistration(i_shackDict);
		ERROR_FORWARD(*tsm);

		///sweet::Shacks::ProgramArgumentsDictionary* shackDict = (sweet::Shacks::ProgramArgumentsDictionary*) i_shackDict;
		bool retval = tsm->setup_3_timestepper(
						i_timestepping_method,
						config
					);

		if (!retval)
		{
			SWEETErrorFatal("Invalid time stepping method: " + i_timestepping_method);
			////int helpVerbosity = 0;
			////if (shackTimeDisc->timestepping_method == "helpall")
			////	helpVerbosity = 1;

			////timeSteppersNewTS.outputHelp(std::cout, "", helpVerbosity);

			////error.clear();
			////std::cout << "Finishing now..." << std::endl;
		}

		tsm->timeIntegrator->setTimeStepSize(i_dt);
		ERROR_FORWARD(*tsm);
	}

public:

	virtual
	bool runIntegration(
				const sweet::XBraid::Vector* i_U,
				const sweet::XBraid::Vector* o_U,
				double i_simulationTime
	)
	{

		ODE_Generic::XBraid::DataContainer* U1 = (ODE_Generic::XBraid::DataContainer*) i_U;
		ODE_Generic::XBraid::DataContainer* U2 = (ODE_Generic::XBraid::DataContainer*) o_U;

		tsm->runIntegration(
				*U1->data,
				*U2->data,
				i_simulationTime
			);
		///U_level->data_DE->swap(*U_level->data_DE_tmp);

		return true;
	}

public:
	bool storePreviousSolution(
				const sweet::XBraid::Vector* i_U
	) override
	{
		ODE_Generic::XBraid::DataContainer* U = (ODE_Generic::XBraid::DataContainer*) i_U;

		tsm->timeIntegrator->storePrevSolution(U->data);

		return true;
	}

};

}}

#endif
