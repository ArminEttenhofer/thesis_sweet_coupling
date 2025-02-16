/*
 * 		Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONCART2D_PROGRAMPDEADVECTIONCART2D_HPP
#define PROGRAMS_PDE_ADVECTIONCART2D_PROGRAMPDEADVECTIONCART2D_HPP


// This is just for the editor to show code as used within precompiler #if ... directives
#include <programs/PDE_AdvectionCart2D/PDEAdvectionCart2DBenchmarksCombined.hpp>
#include <programs/PDE_AdvectionCart2D/PDEAdvectionCart2DTimeSteppers.hpp>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>

// Benchmarks

// Time steppers

#if SWEET_GUI
	#include <sweet/GUI/VisSweet.hpp>
#endif


class ProgramPDEAdvectionCart2D
#if SWEET_GUI
		:	public sweet::GUI::SimulationGUICallbacks
#endif
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

		sweet::Data::Cart2D::DataSpectral prog_h;
		sweet::Data::Cart2D::DataSpectral prog_h_t0;	// at t0
		sweet::Data::Cart2D::DataSpectral prog_u;
		sweet::Data::Cart2D::DataSpectral prog_v;


		bool setup(sweet::Data::Cart2D::Shack *i_shackCart2DDataOps)
		{
			/*
			 * Setup Cart2D Data Config & Operators
			 */
			cart2DDataConfig.setupAuto(*i_shackCart2DDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(cart2DDataConfig);

			ops.setup(cart2DDataConfig, *i_shackCart2DDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(ops);

			prog_h.setup(cart2DDataConfig);
			prog_h_t0.setup(cart2DDataConfig);
			prog_u.setup(cart2DDataConfig);
			prog_v.setup(cart2DDataConfig);

			return true;
		}

		void clear()
		{
			prog_h.clear();
			prog_h_t0.clear();
			prog_u.clear();
			prog_v.clear();

			ops.clear();
			cart2DDataConfig.clear();
		}
	};

	// Simulation data
	DataConfigOps dataConfigOps;

	// time integrators
	PDEAdvectionCart2DTimeSteppers timeSteppers;

	// Handler to all benchmarks
	PDEAdvectionCart2DBenchmarksCombined cart2dBenchmarksCombined;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict;
	sweet::Data::Cart2D::Shack *shackCart2DDataOps;
	sweet::IO::Shack *shackIOData;
	sweet::TimeTree::Shack *shackTimestepControl;
	ShackPDEAdvectionCart2DTimeDiscretization *shackTimeDisc;


#if SWEET_GUI
	// Data to visualize is stored to this variable
	sweet::Data::Cart2D::DataGrid vis_cart2d_data;

	// Which primitive to use for rendering
	int vis_render_type_of_primitive_id = 0;

	// Which primitive to use for rendering
	int vis_data_id = 0;

#endif

public:
	ProgramPDEAdvectionCart2D(
			int i_argc,
			char *const * const i_argv
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackCart2DDataOps(nullptr),
		shackIOData(nullptr),
		shackTimestepControl(nullptr),
		shackTimeDisc(nullptr)
	{
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN(shackProgArgDict);
	}



	bool setup_1_shackRegistration()
	{
		/*
		 * SHACK: Register classes which we require
		 */
		shackCart2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Cart2D::Shack>();
		shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::TimeTree::Shack>();
		shackIOData = shackProgArgDict.getAutoRegistration<sweet::IO::Shack>();
		shackTimeDisc = shackProgArgDict.getAutoRegistration<ShackPDEAdvectionCart2DTimeDiscretization>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * SHACK: Register other things before parsing program arguments
		 */
		cart2dBenchmarksCombined.shackRegistration(shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(cart2dBenchmarksCombined);

		timeSteppers.shackRegistration(shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);

		shackProgArgDict.processHelpArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		return true;
	}

	void clear_1_shackRegistration()
	{
		shackCart2DDataOps = nullptr;
		shackTimestepControl = nullptr;
		shackIOData = nullptr;
		shackTimeDisc = nullptr;

		cart2dBenchmarksCombined.clear();
		timeSteppers.clear();
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

	bool setup_3_dataOpsEtc()
	{
		/*
		 * Setup Cart2D Data Config & Operators
		 */
		dataConfigOps.setup(shackCart2DDataOps);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataConfigOps);

		/*
		 * After we setup the cart2d, we can setup the time steppers and their buffers
		 */
		timeSteppers.setup(shackProgArgDict, dataConfigOps.ops);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);

#if SWEET_GUI
		vis_cart2d_data.setup(dataConfigOps.cart2DDataConfig);
#endif

		std::cout << "Printing shack information:" << std::endl;
		shackProgArgDict.printShackData();

		cart2dBenchmarksCombined.setupInitialConditions(
				dataConfigOps.prog_h,
				dataConfigOps.prog_u,
				dataConfigOps.prog_v,
				dataConfigOps.ops,
				shackProgArgDict
			);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(cart2dBenchmarksCombined);

		dataConfigOps.prog_h_t0 = dataConfigOps.prog_h;

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
#if SWEET_GUI
		vis_cart2d_data.clear();
#endif

		timeSteppers.clear();

		dataConfigOps.clear();
	}

	bool setup()
	{
		if (!setup_1_shackRegistration())
			return false;

		if (!setup_2_processArguments())
			return false;

		if (!setup_3_dataOpsEtc())
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
		std::cout << "Lmax error: " << (dataConfigOps.prog_h_t0-dataConfigOps.prog_h).toGrid().grid_reduce_max_abs() << std::endl;
		std::cout << "RMS error: " << (dataConfigOps.prog_h_t0-dataConfigOps.prog_h).toGrid().grid_reduce_rms() << std::endl;
	}

	double getErrorLMaxOnH()
	{
		return (dataConfigOps.prog_h_t0-dataConfigOps.prog_h).toGrid().grid_reduce_max_abs();
	}

	virtual ~ProgramPDEAdvectionCart2D()
	{
		clear();
	}

	bool runTimestep()
	{
		shackTimestepControl->timestepHelperStart();

		timeSteppers.master->runTimestep(
				dataConfigOps.prog_h, dataConfigOps.prog_u, dataConfigOps.prog_v,
				shackTimestepControl->currentTimestepSize,
				shackTimestepControl->currentSimulationTime
			);

		shackTimestepControl->timestepHelperEnd();

		if (shackIOData->verbosity > 10)
		{
			std::cout << "ts_nr=" << shackTimestepControl->currentTimestepNr << ", t=" << shackTimestepControl->currentSimulationTime*shackIOData->outputFormatTimeScaleInv << std::endl;
			std::cout << "error:" << getErrorLMaxOnH() << std::endl;
		}
		return true;
	}


	bool should_quit()
	{
		return shackTimestepControl->isFinalTimestepReached();
	}


#if SWEET_GUI
	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(
		int i_num_iterations
	)
	{
		if (shackTimestepControl->run_simulation_timesteps)
			for (int i = 0; i < i_num_iterations && !should_quit(); i++)
				runTimestep();
	}

	void vis_getDataArray(
			const sweet::Data::Cart2D::DataGrid **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive_id,
			void **o_bogus_data,
			double *o_vis_min,
			double *o_vis_max,
			bool *vis_reset
	)
	{
		*o_render_primitive_id = vis_render_type_of_primitive_id;

		int id = vis_data_id % 3;
		switch (id)
		{
		case 0:
			sweet::Data::Convert::Cart2DDataSpectral_2_Cart2DDataGrid::convert(dataConfigOps.prog_h, vis_cart2d_data);
			break;

		case 1:
			sweet::Data::Convert::Cart2DDataSpectral_2_Cart2DDataGrid::convert(dataConfigOps.prog_u, vis_cart2d_data);
			break;

		case 2:
			sweet::Data::Convert::Cart2DDataSpectral_2_Cart2DDataGrid::convert(dataConfigOps.prog_v, vis_cart2d_data);
			break;
		}

		*o_dataArray = &vis_cart2d_data;
		*o_aspect_ratio = 1;
	}


	const std::string vis_getStatusString(bool &o_replace_commas_with_newline)
	{
		const char* description = "";
		int id = vis_data_id % 3;

		switch (id)
		{
		default:
		case 0:
			description = "H";
			break;

		case 1:
			description = "u";
			break;

		case 2:
			description = "v";
			break;
		}

		static char title_string[2048];

		sprintf(title_string,
#if SWEET_MPI
				"Rank %i - "
#endif
				"Time: %f (%.2f d), k: %i, dt: %.3e, Vis: %s, MaxVal: %.6e, MinVal: %.6e ",
#if SWEET_MPI
				-1,	// TODO: mpi_rank,
#endif
				shackTimestepControl->currentSimulationTime,
				shackTimestepControl->currentSimulationTime/(60.0*60.0*24.0),
				shackTimestepControl->currentTimestepNr,
				shackTimestepControl->currentTimestepSize,
				description,
				vis_cart2d_data.grid_reduce_max(),
				vis_cart2d_data.grid_reduce_min()
		);

		return title_string;
	}


	void vis_pause()
	{
		shackTimestepControl->run_simulation_timesteps = !shackTimestepControl->run_simulation_timesteps;
	}


	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			vis_data_id++;
			break;

		case 'V':
			vis_data_id--;
			break;

		case 'b':
			vis_render_type_of_primitive_id = (vis_render_type_of_primitive_id + 1) % 2;
			break;
		}
	}
#endif
};




#endif
