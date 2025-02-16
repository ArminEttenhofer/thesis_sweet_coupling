/*
 * 		Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef PROGRAMS_ODE_SCALAR_PROGRAMXBRAIDODESCALAR_HPP
#define PROGRAMS_ODE_SCALAR_PROGRAMXBRAIDODESCALAR_HPP


// This is just for the editor to show code as used within precompiler #if ... directives
#include <programs/ODE_Scalar/ODEScalarBenchmarksCombined.hpp>
#include <programs/ODE_Scalar/ODEScalarTimeSteppers.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>
#include <sweet/XBraid/Shack.hpp>
#include <sweet/XBraid/XBraid_sweet_lib.hpp>


class ProgramXBraidODEScalar
{
public:
	sweet::Error::Base error;

	/*
	 * Just a class to store simulation data all together
	 */
	class Data
	{
	public:
		sweet::Error::Base error;

		double prog_u;
		double prog_u_t0;

		bool setup()
		{
			return true;
		}

		void clear()
		{
		}
	};

	// Simulation data
	Data data;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict;
	sweet::IO::Shack *shackIOData;
	sweet::TimeTree::Shack *shackTimestepControl;
	ShackODEScalarTimeDiscretization *shackTimeDisc;
	sweet::Parallelization::Shack *shackParallelization;
	ShackODEScalar *shackODEScalar;
	ShackODEScalarBenchmarks *shackBenchmarks;
	sweet::XBraid::Shack *shackXBraid;

	// XBraid
	sweet::XBraid::sweet_BraidApp* xbraid_app = nullptr;
	BraidCore* xbraid_core = nullptr;

	// MPI
	MPI_Comm mpi_comm;
	int mpi_rank;

public:
	ProgramXBraidODEScalar(
			int i_argc,
			char *const * const i_argv,
			MPI_Comm i_mpi_comm,
			int i_mpi_rank
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackIOData(nullptr),
		shackTimestepControl(nullptr),
		shackTimeDisc(nullptr),
		shackParallelization(nullptr),
		shackODEScalar(nullptr),
		shackBenchmarks(nullptr),
		shackXBraid(nullptr),
		mpi_comm(i_mpi_comm),
		mpi_rank(i_mpi_rank)
	{
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN(shackProgArgDict);
	}

	bool setup_1_shackRegistration()
	{
		/*
		 * SHACK: Register classes which we require
		 */
		shackIOData = shackProgArgDict.getAutoRegistration<sweet::IO::Shack>();
		shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::TimeTree::Shack>();
		shackTimeDisc = shackProgArgDict.getAutoRegistration<ShackODEScalarTimeDiscretization>();
		shackParallelization = shackProgArgDict.getAutoRegistration<sweet::Parallelization::Shack>();
		shackODEScalar = shackProgArgDict.getAutoRegistration<ShackODEScalar>();
		shackBenchmarks = shackProgArgDict.getAutoRegistration<ShackODEScalarBenchmarks>();
		shackXBraid = shackProgArgDict.getAutoRegistration<sweet::XBraid::Shack>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackProgArgDict.processHelpArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		return true;
	}

	void clear_1_shackRegistration()
	{
		shackIOData = nullptr;
		shackTimestepControl = nullptr;
		shackTimeDisc = nullptr;
		shackParallelization = nullptr;
		shackODEScalar = nullptr;
		shackBenchmarks = nullptr;
		shackXBraid = nullptr;

		shackProgArgDict.clear();
	}

	bool setup_2_processArguments()
	{
		shackProgArgDict.setup();

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

	bool setup_3_data()
	{

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
		xbraid_app = new sweet::XBraid::sweet_BraidApp(mpi_comm, mpi_rank, 0., shackTimestepControl->maxSimulationTime, nt);//, &shackProgArgDict);
		xbraid_app->shackRegistration(shackProgArgDict);

		// XBraid core
		if (shackXBraid->xbraid_run_wrapper_tests)
			xbraid_app->setup();
		else
		{
			xbraid_core = new BraidCore(mpi_comm, xbraid_app);
			xbraid_app->setup(*xbraid_core);
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

		data.clear();
	}

	bool setup()
	{
		if (!setup_1_shackRegistration())
			return false;

		if (!setup_2_processArguments())
			return false;

		if (!setup_3_data())
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
		std::cout << "Error compared to initial condition" << std::endl;
		std::cout << "Error: " << std::abs(data.prog_u_t0-data.prog_u) << std::endl;
	}

	~ProgramXBraidODEScalar()
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




#endif
