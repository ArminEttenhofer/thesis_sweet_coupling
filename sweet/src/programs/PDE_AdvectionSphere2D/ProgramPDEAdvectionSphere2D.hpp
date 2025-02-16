/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONSPHERE2D_PROGRAMPDEADVECTIONSPHERE2D_HPP
#define PROGRAMS_PDE_ADVECTIONSPHERE2D_PROGRAMPDEADVECTIONSPHERE2D_HPP


// This is just for the editor to show code as used within precompiler #if ... directives
#include <programs/PDE_AdvectionSphere2D/benchmarks/ShackPDEAdvectionSphere2DBenchmarks.hpp>
#include <programs/PDE_AdvectionSphere2D/PDEAdvectionSphere2DBenchmarksCombined.hpp>
#include <programs/PDE_AdvectionSphere2D/PDEAdvectionSphere2DTimeSteppers.hpp>
#include <programs/PDE_AdvectionSphere2D/ShackPDEAdvectionSphere2D.hpp>
#include <sweet/Data/Sphere2D/Shack.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>


#if SWEET_GUI
	#include <sweet/GUI/VisSweet.hpp>
	#include <sweet/Data/Cart2D/Cart2D.hpp>
	#include <sweet/Data/Cart2D/Shack.hpp>
#endif

#include <sweet/Data/Sphere2D/Convert/DataSpectral_2_Cart2D_DataGrid.hpp>
#include <sweet/Data/Sphere2D/Convert/DataGrid_2_Cart2D_DataGrid.hpp>



class ProgramPDEAdvectionSphere2D
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

		sweet::Data::Sphere2D::Config sphere2DDataConfig;
		sweet::Data::Sphere2D::Operators ops;

		std::vector<sweet::Data::Sphere2D::DataSpectral> prog_vec;
		sweet::Data::Sphere2D::DataGrid vel_u;
		sweet::Data::Sphere2D::DataGrid vel_v;

		std::vector<sweet::Data::Sphere2D::DataSpectral> prog_vec_t0;

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
				sweet::Data::Sphere2D::Shack *i_shackSphere2DDataOps,
				int i_num_vec_elements
		)
		{
			/*
			 * Setup Sphere2D Data Config & Operators
			 */
			sphere2DDataConfig.setupAuto(i_shackSphere2DDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphere2DDataConfig);

			ops.setup(&sphere2DDataConfig, i_shackSphere2DDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(ops);

			vel_u.setup(sphere2DDataConfig);
			vel_v.setup(sphere2DDataConfig);

			prog_vec.resize(i_num_vec_elements);
			prog_vec_t0.resize(i_num_vec_elements);
			for (int i = 0; i < i_num_vec_elements; i++)
			{
				prog_vec[i].setup(sphere2DDataConfig);
				prog_vec_t0[i].setup(sphere2DDataConfig);
			}

#if SWEET_GUI
			sweet::Data::Cart2D::Shack shackCart2DDataOps;
			shackCart2DDataOps.space_res_physical[0] = i_shackSphere2DDataOps->space_res_physical[0];
			shackCart2DDataOps.space_res_physical[1] = i_shackSphere2DDataOps->space_res_physical[1];
			shackCart2DDataOps.reuse_spectral_transformation_plans = i_shackSphere2DDataOps->reuse_spectral_transformation_plans;

			cart2DDataConfig.setupAuto(shackCart2DDataOps);
#endif
			return true;
		}

		void clear()
		{
			for (std::size_t i = 0; i < prog_vec.size(); i++)
			{
				prog_vec[i].clear();
				prog_vec_t0[i].clear();
			}

			prog_vec.clear();
			prog_vec_t0.clear();

			vel_u.clear();
			vel_v.clear();

			ops.clear();
			sphere2DDataConfig.clear();
		}
	};

	// Simulation data
	DataConfigOps dataConfigOps;

	// time integrators
	PDEAdvectionSphere2DTimeSteppers timeSteppers;

	// Handler to all benchmarks
	PDEAdvectionSphere2DBenchmarksCombined benchmarksCombined;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict;
	sweet::Data::Sphere2D::Shack *shackSphere2DDataOps;
	sweet::IO::Shack *shackIOData;
	sweet::TimeTree::Shack *shackTimestepControl;
	ShackPDEAdvectionSphere2DTimeDiscretization *shackTimeDisc;
	ShackPDEAdvectionSphere2DBenchmarks *shackBenchmarks;
	ShackPDEAdvectionSphere2D *shackPDEAdvectionSphere2D;



public:
	ProgramPDEAdvectionSphere2D(
			int i_argc,
			char *const * const i_argv
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackSphere2DDataOps(nullptr),
		shackIOData(nullptr),
		shackTimestepControl(nullptr),
		shackTimeDisc(nullptr),
		shackBenchmarks(nullptr),
		shackPDEAdvectionSphere2D(nullptr)
	{
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN(shackProgArgDict);
	}


	bool setup_1_shackRegistration()
	{
		/*
		 * Setup argument parsing
		 */
		shackProgArgDict.setup();

		/*
		 * SHACK: Register classes which we require
		 */
		shackSphere2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Sphere2D::Shack>();
		shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::TimeTree::Shack>();
		shackIOData = shackProgArgDict.getAutoRegistration<sweet::IO::Shack>();
		shackTimeDisc = shackProgArgDict.getAutoRegistration<ShackPDEAdvectionSphere2DTimeDiscretization>();
		shackBenchmarks = shackProgArgDict.getAutoRegistration<ShackPDEAdvectionSphere2DBenchmarks>();
		shackPDEAdvectionSphere2D = shackProgArgDict.getAutoRegistration<ShackPDEAdvectionSphere2D>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * SHACK: Register other things before parsing program arguments
		 */

		/*
		 * Setup benchmarks
		 */
		benchmarksCombined.setup_1_registerAllBenchmark();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(benchmarksCombined);

		benchmarksCombined.setup_2_shackRegistration(&shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(benchmarksCombined);


		/*
		 * Setup time steppers
		 */
		timeSteppers.setup_1_registerAllTimesteppers();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);

		timeSteppers.setup_2_shackRegistration(&shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);

		/*
		 * Process HELP arguments
		 */
		shackProgArgDict.processHelpArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * Close shack registration & getting shacks
		 */
		shackProgArgDict.closeRegistration();
		shackProgArgDict.closeGet();

		return true;
	}

	void clear_1_shackRegistration()
	{
		shackSphere2DDataOps = nullptr;
		shackTimestepControl = nullptr;
		shackIOData = nullptr;
		shackTimeDisc = nullptr;

		benchmarksCombined.clear();
		timeSteppers.clear();
		shackProgArgDict.clear();
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

	void clear_2_processArguments()
	{
		shackProgArgDict.clear();
	}

	bool setup_3_dataOpsEtc()
	{
		/*
		 * BENCHMARK: Detect particular benchmark to use
		 */
		benchmarksCombined.setup_3_benchmarkDetection();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(benchmarksCombined);

		/*
		 * Setup benchmark itself
		 */
		benchmarksCombined.setup_4_benchmarkSetup_1_withoutOps();

		/*
		 * Setup the data fields
		 */
		dataConfigOps.setup(shackSphere2DDataOps, benchmarksCombined.benchmark->getNumPrognosticFields());
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataConfigOps);

		/*
		 * Setup benchmark itself
		 */
		benchmarksCombined.setup_5_benchmarkSetup_2_withOps(&dataConfigOps.ops);

		/*
		 * Now we're ready to setup the time steppers
		 */
		timeSteppers.setup_3_timestepper(
				shackTimeDisc->timestepping_method,
				&shackProgArgDict,
				&dataConfigOps.ops
			);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);


		/*
		 * Load initial state of benchmark
		 */
		benchmarksCombined.benchmark->getInitialState(
				dataConfigOps.prog_vec,
				dataConfigOps.vel_u,
				dataConfigOps.vel_v
			);

		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(benchmarksCombined);


		/*
		 * Backup data at t=0
		 */
		dataConfigOps.prog_vec_t0 = dataConfigOps.prog_vec;

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
		dataConfigOps.vis_cart2d_data.clear();
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

		std::cout << "Printing shack information:" << std::endl;
		shackProgArgDict.printShackData();

		std::cout << "SETUP FINISHED" << std::endl;
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
		clear();

		if (!setup())
		{
			error.print();
			return false;
		}

		return !error.exists();
	}

	double getErrorLMaxOnH()
	{
		return (dataConfigOps.prog_vec_t0[0]-dataConfigOps.prog_vec[0]).toGrid().grid_reduce_max_abs();
	}

	double getErrorRMSOnH()
	{
		return (dataConfigOps.prog_vec_t0[0]-dataConfigOps.prog_vec[0]).toGrid().grid_reduce_rms();
	}

	void printSimulationErrors()
	{
		std::cout << "Error compared to initial condition" << std::endl;
		std::cout << "Lmax error: " << getErrorLMaxOnH() << std::endl;
		std::cout << "RMS error: " << getErrorRMSOnH() << std::endl;
	}

	virtual ~ProgramPDEAdvectionSphere2D()
	{
		clear();
	}


	bool runTimestep()
	{
		shackTimestepControl->timestepHelperStart();


		timeSteppers.timestepper->runTimestep(
				dataConfigOps.prog_vec, dataConfigOps.vel_u, dataConfigOps.vel_v,
				shackTimestepControl->currentTimestepSize,
				shackTimestepControl->currentSimulationTime
			);

		if (shackBenchmarks->getVelocities)
		{
			/*
			 * Update velocities just for sake of the correction visualization
			 */
			shackBenchmarks->getVelocities(
					dataConfigOps.vel_u,
					dataConfigOps.vel_v,
					shackTimestepControl->currentSimulationTime + shackTimestepControl->currentTimestepSize,
					shackBenchmarks->getVelocitiesUserData
				);
		}

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
	void vis_post_frame_processing(int i_num_iterations)
	{
		if (shackTimestepControl->run_simulation_timesteps)
			for (int i = 0; i < i_num_iterations && !should_quit(); i++)
				runTimestep();
	}


	void vis_getDataArray(
			const sweet::Data::Cart2D::DataGrid **o_dataArray,
			double *o_aspect_ratio,
			int *o_vis_render_type_of_primitive_id,
			void **o_bogus_data,
			double *o_viz_min,
			double *o_viz_max,
			bool *viz_reset
	)
	{
		*o_vis_render_type_of_primitive_id = dataConfigOps.vis_render_type_of_primitive_id;
		*o_bogus_data = &dataConfigOps.sphere2DDataConfig;

		if (dataConfigOps.vis_data_id < 0)
		{
			int id = -dataConfigOps.vis_data_id-1;

			if (id <  (int)dataConfigOps.prog_vec.size())
			{
				dataConfigOps.vis_cart2d_data = sweet::Data::Sphere2D::Convert::DataSpectral_2_Cart2D_DataGrid::convert(
							dataConfigOps.prog_vec[id] - dataConfigOps.prog_vec_t0[id],
							dataConfigOps.cart2DDataConfig
						);

				*o_dataArray = &dataConfigOps.vis_cart2d_data;
				*o_aspect_ratio = 0.5;
				return;
			}
		}

		std::size_t id = dataConfigOps.vis_data_id;

		if (id >= 0 && id < dataConfigOps.prog_vec.size())
		{
			dataConfigOps.vis_cart2d_data = sweet::Data::Sphere2D::Convert::DataSpectral_2_Cart2D_DataGrid::convert(dataConfigOps.prog_vec[id], dataConfigOps.cart2DDataConfig);
		}
		else if (id >= dataConfigOps.prog_vec.size() && id < dataConfigOps.prog_vec.size() + 2)
		{
			switch (id - dataConfigOps.prog_vec.size())
			{
			case 0:
				dataConfigOps.vis_cart2d_data = sweet::Data::Convert::Sphere2DDataGrid_2_Cart2DDataGrid::grid_convert(dataConfigOps.vel_u, dataConfigOps.cart2DDataConfig);
				break;

			case 1:
				dataConfigOps.vis_cart2d_data = sweet::Data::Convert::Sphere2DDataGrid_2_Cart2DDataGrid::grid_convert(dataConfigOps.vel_v, dataConfigOps.cart2DDataConfig);
				break;
			}
		}
		else if (
			dataConfigOps.prog_vec.size() == 2		&&
			id >= dataConfigOps.prog_vec.size() + 2	&&
			id < dataConfigOps.prog_vec.size() + 4
		)
		{
			sweet::Data::Sphere2D::DataGrid u, v;
			dataConfigOps.ops.vrtdiv_2_uv(dataConfigOps.prog_vec[0], dataConfigOps.prog_vec[1], u, v);

			switch (id - dataConfigOps.prog_vec.size() - 2)
			{
			case 0:
				dataConfigOps.vis_cart2d_data = sweet::Data::Convert::Sphere2DDataGrid_2_Cart2DDataGrid::grid_convert(u, dataConfigOps.cart2DDataConfig);
				break;

			case 1:
				dataConfigOps.vis_cart2d_data = sweet::Data::Convert::Sphere2DDataGrid_2_Cart2DDataGrid::grid_convert(v, dataConfigOps.cart2DDataConfig);
				break;
			}
		}
		else
		{
			SWEET_ASSERT(dataConfigOps.vis_cart2d_data.grid_space_data != nullptr);
			SWEET_ASSERT(dataConfigOps.vis_cart2d_data.cart2DDataConfig != nullptr);
			dataConfigOps.vis_cart2d_data.grid_setZero();
		}

		*o_dataArray = &dataConfigOps.vis_cart2d_data;
		*o_aspect_ratio = 0.5;
	}



	const std::string vis_getStatusString(bool &o_replace_commas_with_newline)
	{
		std::string description = "";

		bool found = false;
		if (dataConfigOps.vis_data_id < 0)
		{
			int id = -dataConfigOps.vis_data_id-1;

			if (id <  (int)dataConfigOps.prog_vec.size())
			{
				std::ostringstream msg;
				msg << "DIFF prog. field " << id;
				description = msg.str();
			}
		}

		if (!found)
		{
			std::size_t id = dataConfigOps.vis_data_id;

			if (id >= 0 && id < dataConfigOps.prog_vec.size())
			{
				std::ostringstream msg;
				msg << "Prog. field " << id;
				description = msg.str();
			}
			else if (id >= dataConfigOps.prog_vec.size() && id < dataConfigOps.prog_vec.size() + 2)
			{
				switch (id - dataConfigOps.prog_vec.size())
				{
				case 0:
					description = "u velocity";
					break;

				case 1:
					description = "v velocity";
					break;
				}
			}
			else if (
					dataConfigOps.prog_vec.size() == 2		&&
					id >= dataConfigOps.prog_vec.size() + 2	&&
					id < dataConfigOps.prog_vec.size() + 4
			)
			{
				switch (id - dataConfigOps.prog_vec.size() - 2)
				{
				case 0:
					description = "prognostic field: u velocity";
					break;

				case 1:
					description = "prognostic field: v velocity";
					break;
				}
			}
			else
			{
				description = "field doesn't exist";
			}
		}

		static char title_string[2048];

		//sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.14e, Vis: %s, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
		sprintf(title_string,
				"Time: %f (%.2f d), k: %i, dt: %.3e, Vis: %s, MaxVal: %.6e, MinVal: %.6e "
				","
				"Colorscale: lowest [Blue... green ... red] highest",
				shackTimestepControl->currentSimulationTime,
				shackTimestepControl->currentSimulationTime/(60.0*60.0*24.0),
				shackTimestepControl->currentTimestepNr,
				shackTimestepControl->currentTimestepSize,
				description.c_str(),
				dataConfigOps.vis_cart2d_data.grid_reduce_max(),
				dataConfigOps.vis_cart2d_data.grid_reduce_min()
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
			dataConfigOps.vis_data_id++;
			break;

		case 'V':
			dataConfigOps.vis_data_id--;
			break;

		case 'b':
		case 'B':
			dataConfigOps.vis_render_type_of_primitive_id = (dataConfigOps.vis_render_type_of_primitive_id + 1) % 2;
			break;
		}
	}
#endif
};



#endif
