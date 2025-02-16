/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_ODE_GENERIC_PROGRAM_HPP
#define PROGRAMS_ODE_GENERIC_PROGRAM_HPP


#include <sweet/Tools/StopwatchBox.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/Parallelization/Shack.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>


#include "Shack.hpp"
#include "DE_Base.hpp"
#include "DE_Dahlquist/DE_Dahlquist.hpp"


namespace ODE_Generic {

class Program
{
public:
	sweet::Error::Base error;

	std::shared_ptr<ODE_Generic::DE_Base> deBase;

	/*
	 * Just a class to store simulation data all together
	 */
	class DataConfigOps
	{
	public:
		sweet::Error::Base error;

		sweet::Data::GenericContainer::Base *prog;
		sweet::Data::GenericContainer::Base *progTmp;


		DataConfigOps()	:
			prog(nullptr),
			progTmp(nullptr)
		{

		}

		void setup(DE_Base *i_deBase)
		{
			prog = i_deBase->getNewDataContainerInstance();
			progTmp = i_deBase->getNewDataContainerInstance();
		}

		void clear()
		{
			if (prog != nullptr)
			{
				delete prog;
				prog = nullptr;
			}
			if (progTmp != nullptr)
			{
				delete progTmp;
				progTmp = nullptr;
			}
		}

		~DataConfigOps()
		{
			clear();
		}

	};

	// Simulation data
	DataConfigOps dataConfigOps;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict;
	sweet::IO::Shack *shackIOData;
	sweet::TimeTree::Shack *shackTimeTree;
	sweet::Parallelization::Shack *shackParallelization;
	Shack *shackODEGeneric;
	
	//int timestep_nr_last_output_simtime = -1;


public:
	Program(
			int i_argc,
			char *const * const i_argv
	)	:
		deBase(nullptr),
		shackProgArgDict(i_argc, i_argv),
		shackIOData(nullptr),
		shackTimeTree(nullptr),
		shackParallelization(nullptr),
		shackODEGeneric(nullptr)
	{
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN(shackProgArgDict);
	}


	bool setup_1_shackRegistration()
	{
		/*
		 * Setup argument parsing
		 */
		shackProgArgDict.setup();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * SHACK: Register classes which we require
		 */
		shackParallelization = shackProgArgDict.getAutoRegistration<sweet::Parallelization::Shack>();
		shackTimeTree = shackProgArgDict.getAutoRegistration<sweet::TimeTree::Shack>();
		shackIOData = shackProgArgDict.getAutoRegistration<sweet::IO::Shack>();
		shackODEGeneric = shackProgArgDict.getAutoRegistration<Shack>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * Process HELP arguments
		 */
		bool helpRequested = !shackProgArgDict.processHelpArguments();

		if (helpRequested)
		{
			std::cout << "Supported ODEs:" << std::endl;
			std::cout << " + 'dahlquist': Dahlquist equation with 3 linear terms l1, l2, l3" << std::endl;
		}
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		return true;
	}

	void clear_1_shackRegistration()
	{
		shackTimeTree = nullptr;
		shackParallelization = nullptr;
		shackIOData = nullptr;
		shackODEGeneric = nullptr;
	}


	bool setup_2_processArguments()
	{
		/*
		 * SHACK: Process arguments
		 */
		shackProgArgDict.processProgramArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * Do some validation of program arguments
		 */
		shackTimeTree->validateTimestepSize();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*shackTimeTree);

		return true;
	}

	void clear_2_processArguments()
	{
		shackProgArgDict.clear();
	}


	bool setup_3_dataAndOps()
	{
		/*
		 * SHACK: Register other things before parsing program arguments
		 */

		if (shackODEGeneric->ode == "")
			return error.set("Use --ode=... to choose ODE\nSupported ODEs: 'dahlquist'");

		/*
		 * That's the point where you can add other ODEs
		 */
		if (shackODEGeneric->ode == "dahlquist")
		{
			deBase = std::make_shared<ODE_Generic::DE_Dahlquist::DE_Dahlquist>();
		}
		else
		{
			return error.set("Unknown ODE '"+shackODEGeneric->ode+"'");
		}

		/*
		 * Setup benchmarks
		 */
		deBase->setup(&shackProgArgDict, shackODEGeneric->timestepping_method);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*deBase);

		deBase->setTimeStepSize(shackTimeTree->currentTimestepSize);

		/*
		 * Setup the data fields
		 */
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataConfigOps);
		dataConfigOps.setup(deBase.get());


		shackProgArgDict.closeRegistration();
		shackProgArgDict.closeGet();


		/*
		 * Load initial state of benchmark
		 */
		deBase->getInitialState(
				*(dataConfigOps.prog)
			);

		/*
		 * Finish registration & getting class interfaces so that nobody can do some
		 * strange things with this anymore
		 */
		shackProgArgDict.closeRegistration();
		shackProgArgDict.closeGet();

		/*
		 * Now we should check that all program arguments have really been parsed
		 */
		shackProgArgDict.checkAllArgumentsProcessed();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		return true;
	}

	void clear_3_data()
	{
		dataConfigOps.clear();

		if (deBase)
		{
			deBase->clear();
			deBase.reset();
		}
	}

	bool setup()
	{
		sweet::Tools::StopwatchBox::getInstance().main_setup.start();

		if (!setup_1_shackRegistration())
			return false;

		if (!setup_2_processArguments())
			return false;

		if (!setup_3_dataAndOps())
			return false;

		if (shackParallelization->isMPIRoot)
		{
			std::cout << "Printing shack information:" << std::endl;

			shackProgArgDict.printShackData();

			std::cout << "SETUP FINISHED" << std::endl;
		}

		sweet::Tools::StopwatchBox::getInstance().main_setup.stop();

		return true;
	}

	void clear()
	{
		clear_3_data();
		clear_2_processArguments();
		clear_1_shackRegistration();
	}

	bool reset()
	{
		// keep pausing simulation
		bool run_simulation_timesteps = shackTimeTree->runSimulationTimesteps;

		clear();

		if (!setup())
		{
			error.print();
			return false;
		}

		shackTimeTree->runSimulationTimesteps = run_simulation_timesteps;

		return !error.exists();
	}

	virtual
	~Program()
	{
		clear();
	}


	bool runTimestep()
	{
		if (shackTimeTree->timestepHelperStart())
		{
			// update time step size if it's changed!
			deBase->setTimeStepSize(shackTimeTree->currentTimestepSize);
		}

		deBase->runTimestep(
				*(dataConfigOps.prog),
				*(dataConfigOps.progTmp),
				shackTimeTree->currentSimulationTime
			);
		dataConfigOps.prog->swap(*(dataConfigOps.progTmp));


		shackTimeTree->timestepHelperEnd();

		if (shackIOData->verbosity > 2)
			if (shackParallelization->isMPIRoot)
			{
				double output_time = shackTimeTree->currentSimulationTime*shackTimeTree->outputSimulationTimeMultiplier;
				std::cout << shackTimeTree->currentTimestepNr << ": " << output_time << std::endl;
			}

		return true;
	}


	/*!
	 * Do some Output to the console or files
	 */
	bool runOutput()
	{
		if (!shackIOData->outputHelper(shackTimeTree))
			return false;

		//deBase->runOutput(shackIOData->);
		return true;
	}


	bool should_quit()
	{
		return shackTimeTree->isFinalTimestepReached();
	}


	void _timestepDoOutput()
	{
		if (shackODEGeneric->compute_errors)
		{
			double error = deBase->computeError(
					*dataConfigOps.prog,
					shackTimeTree->currentSimulationTime,
					"lmax"
				);

			if (shackParallelization->isMPIRoot)
			{
				if (shackTimeTree->currentTimestepNr == 0)
					std::cout << "[MULE] description.error.[timestep_nr]: [time] [error]" << std::endl;

				std::cout << "[MULE] error."
						<< shackTimeTree->currentTimestepNr
						<< ": " << shackTimeTree->currentSimulationTime
						<< "\t" << error << std::endl;
			}
		}

		deBase->runFileOutput(*dataConfigOps.prog);
	}


public:
	bool timestepHandleOutput()
	{
		if (shackIOData->outputHelper(shackTimeTree))
			_timestepDoOutput();

#if 0
		if (shackIOData->outputEachSimTime < 0)
			return false;

		if (shackTimeTree->currentSimulationTime == timestep_nr_last_output_simtime)
			return false;

		timestep_nr_last_output_simtime = shackTimeTree->currentSimulationTime;

		if (shackTimeTree->currentSimulationTime < shackTimeTree->maxSimulationTime - shackIOData->outputEachSimTime*1e-10)
		{
			if (shackIOData->_outputNextSimTime > shackTimeTree->currentSimulationTime)
				return false;
		}

		_timestepDoOutput();

#endif
		return true;
	}
};

}

#endif
