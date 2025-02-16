/*
 * 		Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef PROGRAMS_ODE_GENERIC_PROGRAMXBRAID_HPP
#define PROGRAMS_ODE_GENERIC_PROGRAMXBRAID_HPP

#include <sweet/Data/GenericContainer/CastHelper.hpp>

#include <sweet/Error/Base.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>
#include <sweet/XBraid/Shack.hpp>
#include <sweet/XBraid/App.hpp>

#include <programs/ODE_Generic/Shack.hpp>
#include <programs/ODE_Generic/DE_Base.hpp>
#include <programs/ODE_Generic/DE_Dahlquist/DE_Dahlquist.hpp>

#include <programs/ODE_Generic/XBraid/App.hpp>
#include <programs/ODE_Generic/XBraid/FileOutput.hpp>
#include <programs/ODE_Generic/XBraid/Error.hpp>

#include <programs/ODE_Generic/DE_Dahlquist/XBraid/FileOutput.hpp>

namespace ODE_Generic {

class ProgramXBraid
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

		ODE_Generic::XBraid::DataContainer* data_xbraid;
		////sweet::Data::GenericContainer::Base *data_xbraid;
		////sweet::Data::GenericContainer::Base *progTmp;

		// dummy stuff for XBraid
		//////sweet::Data::Scalar::Config scalarDataConfig;
		//////sweet::Data::Scalar::Operators ops;
		//////sweet::Data::Scalar::Operators opsComplex;

		DataConfigOps()	:
			data_xbraid(nullptr)
			/////progTmp(nullptr)
		{

		}

		bool setup(Shack *i_shackODEGeneric)
		{
			if (i_shackODEGeneric->ode == "dahlquist")
				data_xbraid = new ODE_Generic::DE_Dahlquist::XBraid::DataContainer;
			else
				return error.set("Unknown ODE '"+i_shackODEGeneric->ode+"'");
			///data_xbraid = i_deBase->getNewDataContainerInstance();
			data_xbraid->setup(0);
			////progTmp = i_deBase->getNewDataContainerInstance();

			return true;
		}

		void clear()
		{
			if (data_xbraid != nullptr)
			{
				data_xbraid->clear();
				data_xbraid = nullptr;
			}
			////if (progTmp != nullptr)
			////{
			////	delete progTmp;
			////	progTmp = nullptr;
			////}
		}

		~DataConfigOps()
		{
			clear();
		}

	};

	// Simulation data
	DataConfigOps dataConfigOps;
	DataConfigOps dataConfigOps_initial_guess;

	// XBraid FileOutput and Error
	ODE_Generic::XBraid::FileOutput* file_output = nullptr;
	ODE_Generic::XBraid::Error* error_xbraid = nullptr;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict;
	sweet::IO::Shack *shackIOData;
	sweet::TimeTree::Shack *shackTimestepControl;
	sweet::Parallelization::Shack *shackParallelization;
	sweet::XBraid::Shack *shackXBraid;
	Shack *shackODEGeneric;

	// XBraid
	ODE_Generic::XBraid::App* xbraid_app = nullptr;
	BraidCore* xbraid_core = nullptr;

	// MPI
	MPI_Comm mpi_comm;
	int mpi_rank;

public:
	ProgramXBraid(
			int i_argc,
			char *const * const i_argv,
			MPI_Comm i_mpi_comm,
			int i_mpi_rank
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackIOData(nullptr),
		shackTimestepControl(nullptr),
		shackParallelization(nullptr),
		shackXBraid(nullptr),
		shackODEGeneric(nullptr),
		mpi_comm(i_mpi_comm),
		mpi_rank(i_mpi_rank)
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
		shackIOData = shackProgArgDict.getAutoRegistration<sweet::IO::Shack>();
		shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::TimeTree::Shack>();
		shackParallelization = shackProgArgDict.getAutoRegistration<sweet::Parallelization::Shack>();
		shackXBraid = shackProgArgDict.getAutoRegistration<sweet::XBraid::Shack>();
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
		shackIOData = nullptr;
		shackTimestepControl = nullptr;
		shackParallelization = nullptr;
		shackXBraid = nullptr;
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
		shackTimestepControl->validateTimestepSize();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*shackTimestepControl);

		return true;
	}

	void clear_2_process_arguments()
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

		if (shackODEGeneric->ode == "dahlquist")
			file_output = new ODE_Generic::DE_Dahlquist::XBraid::FileOutput;
		file_output->setup(
					shackIOData,
					shackTimestepControl,
					shackODEGeneric
		);

		error_xbraid = new ODE_Generic::XBraid::Error;
		error_xbraid->setup(
					shackIOData
		);


		/*
		 * Load initial state of benchmark
		 * Two "initial solutions" are required for XBraid:
		 * - Initial solution at t = 0
		 * - Initial guess (random or zero) at t > 0; not used if skip = 0
		 *   for all time steps of the fine discretization
		 */

		/*
		 * That's the point where you can add other ODEs
		 */
		if (shackODEGeneric->ode == "dahlquist")
		{
			deBase = std::make_shared<ODE_Generic::DE_Dahlquist::DE_Dahlquist>();

			deBase->setup(&shackProgArgDict, shackODEGeneric->timestepping_method);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*deBase);
			///deBase->setTimeStepSize(shackTimestepControl->currentTimestepSize);

			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataConfigOps);
			dataConfigOps.setup(shackODEGeneric);
			dataConfigOps_initial_guess.setup(shackODEGeneric);

			if (shackXBraid->xbraid_use_rand)
			{
				//sweet::Data::GenericContainer::CastHelper<ODE_Generic::DE_Dahlquist::DataContainer::Simulation, ODE_Generic::DE_Dahlquist::DataContainer::Config>::cast(*dataConfigOps_initial_guess.data_xbraid).data[0] = ((double)braid_Rand())/braid_RAND_MAX;
				//dataConfigOps_initial_guess.data_xbraid->data[0] = ((double)braid_Rand())/braid_RAND_MAX;
				dataConfigOps_initial_guess.data_xbraid->op_setValue(((double)braid_Rand())/braid_RAND_MAX);
			}
			else
			{
				dataConfigOps_initial_guess.data_xbraid->op_setZero();
			}

		}
		else
		{
			return error.set("Unknown ODE '"+shackODEGeneric->ode+"'");
		}

		// set initial solution to ODE_Generic::DataContainer then to ODE_Generic::XBraid::DataContainer
		sweet::Data::GenericContainer::Base* tmp;
		tmp = deBase->getNewDataContainerInstance();
		deBase->getInitialState(
				*tmp
			);
		if (shackODEGeneric->ode == "dahlquist")
		{
			// ugly!!
			std::complex<double> v = 
				sweet::Data::GenericContainer::CastHelper<ODE_Generic::DE_Dahlquist::DataContainer::Simulation, ODE_Generic::DE_Dahlquist::DataContainer::Config>::cast(*tmp).data[0];
			dataConfigOps.data_xbraid->op_setValue(v);
		}
		delete tmp;

		std::cout << "Printing shack information:" << std::endl;
		shackProgArgDict.printShackData();

		//////////////////
		// SETUP XBRAID //
		//////////////////

		// get the number of timesteps in the finest level
		int nt = (int) (shackTimestepControl->maxSimulationTime / shackTimestepControl->currentTimestepSize);
		if (nt * shackTimestepControl->currentTimestepSize < shackTimestepControl->maxSimulationTime - 1e-10)
			nt++;

		// XBraid app (user-defined)
		xbraid_app = new ODE_Generic::XBraid::App(
							mpi_comm,
							mpi_rank,
							0.,
							shackTimestepControl->maxSimulationTime,
							nt,
							shackODEGeneric->ode
							/////&dataConfigOps.scalarDataConfig,
							/////&dataConfigOps.ops,
							/////&dataConfigOps.opsComplex
		);//, &shackProgArgDict);
		xbraid_app->shackRegistration(shackProgArgDict);

		////// XBraid core
		////if (shackXBraid->xbraid_run_wrapper_tests)
		////	xbraid_app->setup(*dataConfigOps.prog, *dataConfigOps_initial_guess.prog);
		////else
		////{
		////	xbraid_core = new BraidCore(mpi_comm, xbraid_app);
		////	xbraid_app->setup(*xbraid_core, *dataConfigOps.prog, *dataConfigOps_initial_guess.prog);
		////}
		// XBraid core
		if (shackXBraid->xbraid_run_wrapper_tests)
			xbraid_app->setup(*dataConfigOps.data_xbraid, *dataConfigOps_initial_guess.data_xbraid, file_output, error_xbraid);
		else
		{
			xbraid_core = new BraidCore(mpi_comm, xbraid_app);
			xbraid_app->setup(*xbraid_core, *dataConfigOps.data_xbraid, *dataConfigOps_initial_guess.data_xbraid, file_output, error_xbraid);
		}

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

		/////////timeSteppers.clear();

		if (xbraid_core)
		{
			delete xbraid_core;
			xbraid_core = nullptr;
		}

		if (xbraid_app)
		{
			delete xbraid_app;
			xbraid_app = nullptr;
		}

		dataConfigOps.clear();
		dataConfigOps_initial_guess.clear();
	}

	bool setup()
	{
		if (!setup_1_shackRegistration())
			return false;

		if (!setup_2_processArguments())
			return false;

		if (!setup_3_dataAndOps())
			return false;

		std::cout << "SETUP FINISHED" << std::endl;
		return true;
	}
	void clear()
	{
		clear_3_data();
		clear_2_process_arguments();
		clear_1_shackRegistration();
	}

	bool reset()
	{
		clear();

		if (!setup())
		{
			error.print();
			return false;
		}

		return !error.exists();
	}

	void printSimulationErrors()
	{
	}

	~ProgramXBraid()
	{
		clear();
	}


	bool runXBraid()
	{

		shackTimestepControl->timestepHelperStart();

		// Run wrapper tests
		if (shackXBraid->xbraid_run_wrapper_tests)
		{
			BraidUtil braid_util;
			int test = braid_util.TestAll(xbraid_app, mpi_comm, stdout, 0., shackTimestepControl->currentTimestepSize, shackTimestepControl->currentTimestepSize * 2);
			if (test == 0)
				SWEETErrorFatal("Tests failed!");
			else
				std::cout << "Tests successful!" << std::endl;
		}
		else
		{
			// Run Simulation
			xbraid_core->Drive();
		}

		shackTimestepControl->timestepHelperEnd();

		return true;

	}

	bool should_quit()
	{
		return false;
		////////////return shackTimestepControl->isFinalTimestepReached();
	}

};

};




#endif
