/*
 * 		Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_PROGRAMXBRAIDPDESWESPHERE2D_HPP
#define PROGRAMS_PDE_SWESPHERE2D_PROGRAMXBRAIDPDESWESPHERE2D_HPP


// This is just for the editor to show code as used within precompiler #if ... directives
#include <programs/PDE_SWESphere2D/Benchmarks/Shack.hpp>
#include <programs/PDE_SWESphere2D/BenchmarksCombined.hpp>
///#include <programs/PDE_SWESphere2D/TimeSteppers.hpp>
#include <programs/PDE_SWESphere2D/Diagnostics.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Config.hpp>
///#include <sweet/Data/Sphere2D/GridMapping.hpp>
#include <sweet/Data/Sphere2D/Shack.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>
#include <sweet/XBraid/Shack.hpp>
#include <sweet/XBraid/App.hpp>
#include <programs/PDE_SWESphere2D/XBraid/App.hpp>
#include <programs/PDE_SWESphere2D/XBraid/Error.hpp>
#include <programs/PDE_SWESphere2D/XBraid/FileOutput.hpp>
#include <programs/PDE_SWESphere2D/DataContainer/Topography.hpp>

#include<vector>

namespace PDE_SWESphere2D {

namespace XBraid {

class Program
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

		sweet::Data::Sphere2D::Config sphere2DDataConfig;
		sweet::Data::Sphere2D::Operators ops;
		sweet::Data::Sphere2DComplex::Operators opsComplex;

		PDE_SWESphere2D::XBraid::DataContainer data_xbraid;

		PDE_SWESphere2D::DataContainer::Topography topography;

		sweet::Data::Sphere2D::DataSpectral t0_prog_phi_pert;
		sweet::Data::Sphere2D::DataSpectral t0_prog_vrt;
		sweet::Data::Sphere2D::DataSpectral t0_prog_div;

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

			opsComplex.setup(&sphere2DDataConfig, i_shackSphere2DDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(opsComplex);

			data_xbraid.setup(&sphere2DDataConfig, 0);

			topography.setup(&sphere2DDataConfig);

			t0_prog_phi_pert.setup(sphere2DDataConfig);
			t0_prog_vrt.setup(sphere2DDataConfig);
			t0_prog_div.setup(sphere2DDataConfig);

			return true;
		}

		void clear(bool i_clear_spectral_transforms = true)
		{
			data_xbraid.clear();

			t0_prog_phi_pert.clear();
			t0_prog_vrt.clear();
			t0_prog_div.clear();

			topography.clear();

			ops.clear();

			if (i_clear_spectral_transforms)
				sphere2DDataConfig.clear();
		}
	};

	// Simulation data
	DataConfigOps dataConfigOps;
	DataConfigOps dataConfigOps_initial_guess;

	// XBraid FileOutput and Error
	PDE_SWESphere2D::XBraid::FileOutput* file_output = nullptr;
	PDE_SWESphere2D::XBraid::Error* error_xbraid = nullptr;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict;
	sweet::Data::Sphere2D::Shack *shackSphere2DDataOps;
	sweet::IO::Shack *shackIOData;
	sweet::TimeTree::Shack *shackTimestepControl;
	sweet::Parallelization::Shack *shackParallelization;
	PDE_SWESphere2D::Shack *shackPDESWESphere2D;
	PDE_SWESphere2D::Benchmarks::Shack *shackBenchmarks;
	sweet::XBraid::Shack *shackXBraid;

	// Handler to all benchmarks
	PDE_SWESphere2D::Benchmarks::BenchmarksCombined sphere2dBenchmarksCombined;

	// XBraid
	////sweet::XBraid::App* xbraid_app = nullptr;
	PDE_SWESphere2D::XBraid::App* xbraid_app = nullptr;
	BraidCore* xbraid_core = nullptr;

	// MPI
	MPI_Comm mpi_comm;
	int mpi_rank;

public:
	Program(
			int i_argc,
			char *const * const i_argv,
			MPI_Comm i_mpi_comm,
			int i_mpi_rank
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackSphere2DDataOps(nullptr),
		shackIOData(nullptr),
		shackTimestepControl(nullptr),
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
		shackParallelization = shackProgArgDict.getAutoRegistration<sweet::Parallelization::Shack>();
		shackPDESWESphere2D = shackProgArgDict.getAutoRegistration<PDE_SWESphere2D::Shack>();
		shackBenchmarks = shackProgArgDict.getAutoRegistration<PDE_SWESphere2D::Benchmarks::Shack>();
		shackXBraid = shackProgArgDict.getAutoRegistration<sweet::XBraid::Shack>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		sphere2dBenchmarksCombined.setup_1_registerAllBenchmark();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphere2dBenchmarksCombined);

		sphere2dBenchmarksCombined.setup_2_shackRegistration(&shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphere2dBenchmarksCombined);

		shackProgArgDict.processHelpArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		return true;
	}

	void clear_1_shackRegistration()
	{
		shackSphere2DDataOps = nullptr;
		shackIOData = nullptr;
		shackTimestepControl = nullptr;
		shackParallelization = nullptr;
		shackPDESWESphere2D = nullptr;
		shackBenchmarks = nullptr;
		shackXBraid = nullptr;

		sphere2dBenchmarksCombined.clear();

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
		 * BENCHMARK: Detect particular benchmark to use
		 */
		sphere2dBenchmarksCombined.setup_3_benchmarkDetection();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphere2dBenchmarksCombined);

		/*
		 * Setup benchmark itself
		 */
		sphere2dBenchmarksCombined.setup_4_benchmarkSetup_1_withoutOps();

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
		dataConfigOps_initial_guess.setup(shackSphere2DDataOps);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataConfigOps_initial_guess);

		/*
		 * Setup benchmark itself
		 */
		sphere2dBenchmarksCombined.setup_5_benchmarkSetup_2_withOps(&dataConfigOps.ops);

		file_output = new PDE_SWESphere2D::XBraid::FileOutput;
		file_output->setup(
					shackIOData,
					shackTimestepControl,
					shackSphere2DDataOps,
					shackPDESWESphere2D,
					&dataConfigOps.sphere2DDataConfig,
					&dataConfigOps.ops,
					&dataConfigOps.opsComplex
		);

		error_xbraid = new PDE_SWESphere2D::XBraid::Error;
		error_xbraid->setup(
					&dataConfigOps.sphere2DDataConfig,
					shackIOData
		);

		/*
		 * Load topography
		 */
		sphere2dBenchmarksCombined.benchmark->setup_topography();

		/*
		 * Load initial state of benchmark
		 * Two "initial solutions" are required for XBraid:
		 * - Initial solution at t = 0
		 * - Initial guess (random or zero) at t > 0; not used if skip = 0
		 *   for all time steps of the fine discretization
		 */

		sphere2dBenchmarksCombined.benchmark->getInitialState(
				dataConfigOps.data_xbraid.data->phi_pert,
				dataConfigOps.data_xbraid.data->vrt,
				dataConfigOps.data_xbraid.data->div
			);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphere2dBenchmarksCombined);

		// Initial guess
		if (shackXBraid->xbraid_use_rand)
		{

			sweet::Data::Sphere2D::DataGrid t0_prog_phi_pert_phys(dataConfigOps.sphere2DDataConfig);
			sweet::Data::Sphere2D::DataGrid t0_prog_vrt_phys(dataConfigOps.sphere2DDataConfig);
			sweet::Data::Sphere2D::DataGrid t0_prog_div_phys(dataConfigOps.sphere2DDataConfig);

			t0_prog_phi_pert_phys.grid_update_lambda_array(
						[&](int i, int j, double &io_data)
				{
					io_data = (shackPDESWESphere2D->h0 + ((double)braid_Rand())/braid_RAND_MAX) * shackPDESWESphere2D->gravitation;
				}
			);
			t0_prog_vrt_phys.grid_update_lambda_array(
						[&](int i, int j, double &io_data)
				{
					io_data = ((double)braid_Rand())/braid_RAND_MAX * 1e-6;
				}
			);
			t0_prog_div_phys.grid_update_lambda_array(
						[&](int i, int j, double &io_data)
				{
					io_data = ((double)braid_Rand())/braid_RAND_MAX * 1e-6;
				}
			);

			dataConfigOps_initial_guess.data_xbraid.data->phi_pert.loadSphere2DDataGrid(t0_prog_phi_pert_phys);
			dataConfigOps_initial_guess.data_xbraid.data->vrt.loadSphere2DDataGrid(t0_prog_vrt_phys);
			dataConfigOps_initial_guess.data_xbraid.data->div.loadSphere2DDataGrid(t0_prog_div_phys);
		}
		else
			dataConfigOps_initial_guess.data_xbraid.data->op_setZero();

		// get the number of timesteps in the finest level
		int nt = (int) (shackTimestepControl->maxSimulationTime / shackTimestepControl->currentTimestepSize);
		if (nt * shackTimestepControl->currentTimestepSize < shackTimestepControl->maxSimulationTime - 1e-10)
			nt++;

		// XBraid app (user-defined)
		xbraid_app = new PDE_SWESphere2D::XBraid::App(
							mpi_comm,
							mpi_rank,
							0.,
							shackTimestepControl->maxSimulationTime,
							nt,
							&dataConfigOps.sphere2DDataConfig,
							&dataConfigOps.ops,
							&dataConfigOps.opsComplex
		);//, &shackProgArgDict);
		xbraid_app->shackRegistration(shackProgArgDict);

		// XBraid core
		if (shackXBraid->xbraid_run_wrapper_tests)
			xbraid_app->setup(dataConfigOps.data_xbraid, dataConfigOps_initial_guess.data_xbraid, file_output, error_xbraid);
		else
		{
			xbraid_core = new BraidCore(mpi_comm, xbraid_app);
			xbraid_app->setup(*xbraid_core, dataConfigOps.data_xbraid, dataConfigOps_initial_guess.data_xbraid, file_output, error_xbraid);
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

		if (file_output)
		{
			file_output->clear();
			file_output = nullptr;
		}

		if (error_xbraid)
		{
			error_xbraid->clear();
			error_xbraid = nullptr;
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

	~Program()
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

}}

#endif
