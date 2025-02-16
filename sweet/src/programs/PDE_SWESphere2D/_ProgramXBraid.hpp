/*
 * 	Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_PROGRAMXBRAID_HPP
#define PROGRAMS_PDE_SWESPHERE2D_PROGRAMXBRAID_HPP


// This is just for the editor to show code as used within precompiler #if ... directives
#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/Shack.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>
#include <sweet/XBraid/Shack.hpp>
#include <sweet/XBraid/XBraid_sweet_lib.hpp>

#include "Benchmarks/Shack.hpp"

#include "BenchmarksCombined.hpp"

// Time steppers
#include "TimeOld/PDESWESphere2D_TimeSteppers.hpp"

#include <vector>

namespace PDE_SWESphere2D {

class ProgramXBraid
{
public:
	sweet::Error::Base error;

	/*
	 * Just a class to store simulation data all together
	 */
	class DataConfigOps
	{
	public:
		sweet::Error::Base error;

		sweet::Data::Sphere2D::SphereData_Config sphere2DDataConfig;
		sweet::Data::Sphere2D::SphereOperators ops;
		sweet::Data::Sphere2D::SphereOperatorsComplex opsComplex;

		sweet::Data::Sphere2D::DataSpectral prog_phi_pert;
		sweet::Data::Sphere2D::DataSpectral prog_div;
		sweet::Data::Sphere2D::DataSpectral prog_vrt;


		sweet::Data::Sphere2D::DataSpectral t0_prog_phi_pert;
		sweet::Data::Sphere2D::DataSpectral t0_prog_div;
		sweet::Data::Sphere2D::DataSpectral t0_prog_vrt;

		bool setup(
				sweet::Data::Sphere2D::Shack *i_shackSphere2DDataOps,
				bool i_setup_spectral_transforms = true		// for reset()
		)
		{
			/*
			 * Setup Sphere2D Data Config & Operators
			 */
			if (i_setup_spectral_transforms)
			{
				sphere2DDataConfig.setupAuto(i_shackSphere2DDataOps);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphere2DDataConfig);
			}

			ops.setup(&sphere2DDataConfig, i_shackSphere2DDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(ops);

			opsComplex.setup(&sphereDataConfig, i_shackSphere2DDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(opsComplex);

			prog_phi_pert.setup(sphere2DDataConfig);
			prog_div.setup(sphere2DDataConfig);
			prog_vrt.setup(sphere2DDataConfig);

			return true;
		}

		void clear(bool i_clear_spectral_transforms = true)
		{
			prog_phi_pert.clear();
			prog_div.clear();
			prog_vrt.clear();

			t0_prog_phi_pert.clear();
			t0_prog_div.clear();
			t0_prog_vrt.clear();

			ops.clear();

			if (i_clear_spectral_transforms)
				sphere2DDataConfig.clear();
		}
	};


	// Simulation data
	DataConfigOps dataConfigOps;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict;
	sweet::Data::Sphere2D::Shack *shackSphere2DDataOps;
	sweet::IO::Shack *shackIOData;
	sweet::TimeTree::Shack *shackTimestepControl;
	ShackTimeDiscretization *shackTimeDisc;
	sweet::Parallelization::Shack *shackParallelization;
	PDE_SWESphere2D::Shack *shackPDESWESphere2D;
	Benchmarks::Shack *shackBenchmarks;
	sweet::XBraid::Shack *shackXBraid;

	// XBraid
	sweet::XBraid::sweet_BraidApp* xbraid_app = nullptr;
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
		shackSphere2DDataOps(nullptr),
		shackIOData(nullptr),
		shackTimestepControl(nullptr),
		shackTimeDisc(nullptr),
		shackParallelization(nullptr),
		shackPDESWESphere2D(nullptr),
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
		shackSphere2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Sphere2D::Shack>();
		shackIOData = shackProgArgDict.getAutoRegistration<sweet::IO::Shack>();
		shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::TimeTree::Shack>();
		shackTimeDisc = shackProgArgDict.getAutoRegistration<ShackTimeDiscretization>();
		shackParallelization = shackProgArgDict.getAutoRegistration<sweet::Parallelization::Shack>();
		shackPDESWESphere2D = shackProgArgDict.getAutoRegistration<PDE_SWESphere2D::Shack>();
		shackBenchmarks = shackProgArgDict.getAutoRegistration<Benchmarks::Shack>();
		shackXBraid = shackProgArgDict.getAutoRegistration<sweet::XBraid::Shack>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackProgArgDict.processHelpArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		return true;
	}

	void clear_1_shackRegistration()
	{
		shackSphere2DDataOps = nullptr;
		shackIOData = nullptr;
		shackTimestepControl = nullptr;
		shackTimeDisc = nullptr;
		shackParallelization = nullptr;
		shackPDESWESphere2D = nullptr;
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

		/*
		 * Setup Sphere2D Data Config & Operators
		 */
		dataConfigOps.setup(shackSphere2DDataOps);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataConfigOps);

		// get the number of timesteps in the finest level
		int nt = (int) (shackTimestepControl->maxSimulationTime / shackTimestepControl->currentTimestepSize);
		if (nt * shackTimestepControl->currentTimestepSize < shackTimestepControl->maxSimulationTime - 1e-10)
			nt++;

		// XBraid app (user-defined)
		xbraid_app = new sweet::XBraid::sweet_BraidApp(
								mpi_comm,
								mpi_rank,
								0.,
								shackTimestepControl->max_simulation_time,
								nt,
								&dataConfigOps.sphere2DDataConfig,
								&dataConfigOps.ops,
								&dataConfigOps.opsComplex
							);

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

	////void printSimulationErrors()
	////{
	////	std::cout << "Error compared to initial condition" << std::endl;
	////	std::cout << "Error: " << std::abs(data.prog_u_t0-data.prog_u) << std::endl;
	////}

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
		////////////return shackTimestepControl->isFinalTimestepReached();
		return false;
	}

};

}

#endif
