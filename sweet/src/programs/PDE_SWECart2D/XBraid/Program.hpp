/*
 * 		Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef PROGRAMS_PDE_SWECART2D_PROGRAMXBRAIDPDESWECART2D_HPP
#define PROGRAMS_PDE_SWECART2D_PROGRAMXBRAIDPDESWECART2D_HPP


// This is just for the editor to show code as used within precompiler #if ... directives
#include <programs/PDE_SWECart2D/Benchmarks/Shack.hpp>
#include <programs/PDE_SWECart2D/BenchmarksCombined.hpp>
#include <programs/PDE_SWECart2D/TimeSteppers.hpp>
#include <programs/PDE_SWECart2D/Diagnostics.hpp>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Data/Cart2D/Config.hpp>
#include <sweet/Data/Cart2D/GridMapping.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>
#include <sweet/XBraid/Shack.hpp>
#include <sweet/XBraid/App.hpp>
#include <programs/PDE_SWECart2D/XBraid/App.hpp>
#include <programs/PDE_SWECart2D/XBraid/Error.hpp>
#include <programs/PDE_SWECart2D/XBraid/FileOutput.hpp>

#include<vector>

namespace PDE_SWECart2D {

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

		sweet::Data::Cart2D::Config cart2DDataConfig;
		sweet::Data::Cart2D::Operators ops;
		sweet::Data::Cart2DComplex::Operators opsComplex;

		////DataContainer::Simulation prog;
		////DataContainer::Simulation progTmp;

		PDE_SWECart2D::XBraid::DataContainer data_xbraid;

		sweet::Data::Cart2D::DataSpectral t0_prog_h_pert;
		sweet::Data::Cart2D::DataSpectral t0_prog_u;
		sweet::Data::Cart2D::DataSpectral t0_prog_v;

		// Mapping between grids
		sweet::Data::Cart2D::GridMapping gridMapping;

#if SWEET_GUI
		sweet::Data::Cart2D::Config cart2DDataConfig;

		// Data to visualize is stored to this variable
		sweet::Data::Cart2D::DataGrid vis_cart2d_data;

		// Which primitive to use for rendering
		int vis_render_type_of_primitive_id = 1;

		// Which primitive to use for rendering
		int vis_data_id = 0;
#endif

		bool setup(
				sweet::Data::Cart2D::Shack *i_shackCart2DDataOps,
				bool i_setup_spectral_transforms = true		// for reset()
		)
		{
			/*
			 * Setup Sphere2D Data Config & Operators
			 */
			if (i_setup_spectral_transforms)
			{
				cart2DDataConfig.setupAuto(i_shackCart2DDataOps);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(cart2DDataConfig);
			}

			ops.setup(&cart2DDataConfig, i_shackCart2DDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(ops);

			opsComplex.setup(&cart2DDataConfig, i_shackCart2DDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(opsComplex);

			data_xbraid.setup(&cart2DDataConfig, 0);

			t0_prog_h_pert.setup(cart2DDataConfig);
			t0_prog_u.setup(cart2DDataConfig);
			t0_prog_v.setup(cart2DDataConfig);

#if SWEET_GUI
			sweet::Data::Cart2D::Shack shackCart2DDataOps;
#if 0
			// WARNING: We need to use sphere2DDataConfig, since i_shackSphere2DDataOps is not initialized with a reset()
			shackCart2DDataOps.space_res_physical[0] = i_shackSphere2DDataOps->space_res_physical[0];
			shackCart2DDataOps.space_res_physical[1] = i_shackSphere2DDataOps->space_res_physical[1];
#else
			shackCart2DDataOps.space_res_physical[0] = sphere2DDataConfig.grid_num_lon;
			shackCart2DDataOps.space_res_physical[1] = sphere2DDataConfig.grid_num_lat;

#endif
			shackCart2DDataOps.reuse_spectral_transformation_plans = i_shackCart2DDataOps->reuse_spectral_transformation_plans;

			cart2DDataConfig.setupAuto(shackCart2DDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(cart2DDataConfig);
#endif

			if (i_shackCart2DDataOps->space_grid_use_c_staggering)
				gridMapping.setup(i_shackCart2DDataOps, &cart2DDataConfig);

			return true;
		}

		void clear(bool i_clear_spectral_transforms = true)
		{
			data_xbraid.clear();

			t0_prog_h_pert.clear();
			t0_prog_u.clear();
			t0_prog_v.clear();

			ops.clear();

			if (i_clear_spectral_transforms)
				cart2DDataConfig.clear();
		}
	};

	// Simulation data
	DataConfigOps dataConfigOps;
	DataConfigOps dataConfigOps_initial_guess;

	// XBraid FileOutput and Error
	PDE_SWECart2D::XBraid::FileOutput* file_output = nullptr;
	PDE_SWECart2D::XBraid::Error* error_xbraid = nullptr;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict;
	sweet::Data::Cart2D::Shack *shackCart2DDataOps;
	sweet::IO::Shack *shackIOData;
	sweet::TimeTree::Shack *shackTimestepControl;
	PDE_SWECart2D::TimeDiscretization::Shack *shackTimeDisc;
	sweet::Parallelization::Shack *shackParallelization;
	PDE_SWECart2D::Shack *shackPDESWECart2D;
	PDE_SWECart2D::Benchmarks::Shack *shackBenchmarks;
	sweet::XBraid::Shack *shackXBraid;

	// Handler to all benchmarks
	PDE_SWECart2D::Benchmarks::BenchmarksCombined cart2dBenchmarksCombined;

	// XBraid
	////sweet::XBraid::App* xbraid_app = nullptr;
	PDE_SWECart2D::XBraid::App* xbraid_app = nullptr;
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
		shackCart2DDataOps(nullptr),
		shackIOData(nullptr),
		shackTimestepControl(nullptr),
		shackTimeDisc(nullptr),
		shackParallelization(nullptr),
		shackPDESWECart2D(nullptr),
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
		shackCart2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Cart2D::Shack>();
		shackIOData = shackProgArgDict.getAutoRegistration<sweet::IO::Shack>();
		shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::TimeTree::Shack>();
		shackTimeDisc = shackProgArgDict.getAutoRegistration<PDE_SWECart2D::TimeDiscretization::Shack>();
		shackParallelization = shackProgArgDict.getAutoRegistration<sweet::Parallelization::Shack>();
		shackPDESWECart2D = shackProgArgDict.getAutoRegistration<PDE_SWECart2D::Shack>();
		shackBenchmarks = shackProgArgDict.getAutoRegistration<PDE_SWECart2D::Benchmarks::Shack>();
		shackXBraid = shackProgArgDict.getAutoRegistration<sweet::XBraid::Shack>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		cart2dBenchmarksCombined.shackRegistration(shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(cart2dBenchmarksCombined);

		shackProgArgDict.processHelpArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		return true;
	}

	void clear_1_shackRegistration()
	{
		shackCart2DDataOps = nullptr;
		shackIOData = nullptr;
		shackTimestepControl = nullptr;
		shackTimeDisc = nullptr;
		shackParallelization = nullptr;
		shackPDESWECart2D = nullptr;
		shackBenchmarks = nullptr;
		shackXBraid = nullptr;

		cart2dBenchmarksCombined.clear();

		/////scalarBenchmarksCombined.clear();
		/////timeSteppers.clear();
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
		 * Setup Cart2D Data Config & Operators
		 */
		dataConfigOps.setup(shackCart2DDataOps);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataConfigOps);
		dataConfigOps_initial_guess.setup(shackCart2DDataOps);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataConfigOps_initial_guess);

		file_output = new PDE_SWECart2D::XBraid::FileOutput;
		file_output->setup(
					shackIOData,
					shackTimestepControl,
					shackCart2DDataOps,
					shackPDESWECart2D,
					&dataConfigOps.cart2DDataConfig,
					&dataConfigOps.ops,
					&dataConfigOps.opsComplex
		);

		error_xbraid = new PDE_SWECart2D::XBraid::Error;
		error_xbraid->setup(
					&dataConfigOps.cart2DDataConfig,
					shackIOData
		);

		/*
		 * Load initial state of benchmark
		 * Two "initial solutions" are required for XBraid:
		 * - Initial solution at t = 0
		 * - Initial guess (random or zero) at t > 0; not used if skip = 0
		 *   for all time steps of the fine discretization
		 */

		cart2dBenchmarksCombined.setupInitialConditions(
				dataConfigOps.data_xbraid.data->h_pert,
				dataConfigOps.data_xbraid.data->u,
				dataConfigOps.data_xbraid.data->v,
				&dataConfigOps.ops,
				&dataConfigOps.cart2DDataConfig
			);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(cart2dBenchmarksCombined);

		// Initial guess
		if (shackXBraid->xbraid_use_rand)
		{

			sweet::Data::Cart2D::DataGrid t0_prog_h_phys(&dataConfigOps.cart2DDataConfig);
			sweet::Data::Cart2D::DataGrid t0_prog_u_phys(&dataConfigOps.cart2DDataConfig);
			sweet::Data::Cart2D::DataGrid t0_prog_v_phys(&dataConfigOps.cart2DDataConfig);

			t0_prog_h_phys.grid_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
				{
					io_data = shackPDESWECart2D->h0 + ((double)braid_Rand())/braid_RAND_MAX;
				}
			);
			t0_prog_u_phys.grid_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
				{
					io_data = ((double)braid_Rand())/braid_RAND_MAX;
				}
			);
			t0_prog_v_phys.grid_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
				{
					io_data = ((double)braid_Rand())/braid_RAND_MAX;
				}
			);

			dataConfigOps_initial_guess.data_xbraid.data->h_pert.loadCart2DDataGrid(t0_prog_h_phys);
			dataConfigOps_initial_guess.data_xbraid.data->u.loadCart2DDataGrid(t0_prog_u_phys);
			dataConfigOps_initial_guess.data_xbraid.data->v.loadCart2DDataGrid(t0_prog_v_phys);
		}
		else
			dataConfigOps_initial_guess.data_xbraid.data->op_setZero();


		// get the number of timesteps in the finest level
		int nt = (int) (shackTimestepControl->maxSimulationTime / shackTimestepControl->currentTimestepSize);
		if (nt * shackTimestepControl->currentTimestepSize < shackTimestepControl->maxSimulationTime - 1e-10)
			nt++;

		// XBraid app (user-defined)
		xbraid_app = new PDE_SWECart2D::XBraid::App(
							mpi_comm,
							mpi_rank,
							0.,
							shackTimestepControl->maxSimulationTime,
							nt,
							&dataConfigOps.cart2DDataConfig,
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
