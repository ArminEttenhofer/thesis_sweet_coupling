/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_SCONS_OPTIONS: --sphere2d-spectral-space=enable
 * MULE_SCONS_OPTIONS: --gui=enable
 */


// This is just for the editor to show code as used within precompiler #if ... directives
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Data/Sphere2D/Shack.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>

#if SWEET_GUI
	#include <sweet/GUI/VisSweet.hpp>
	#include <sweet/Data/Cart2D/Cart2D.hpp>
	#include <sweet/Data/Cart2D/Shack.hpp>
#endif

#include <sweet/Data/Sphere2D/Convert/DataSpectral_2_Cart2D_DataGrid.hpp>
#include <sweet/Data/Sphere2D/Convert/DataGrid_2_Cart2D_DataGrid.hpp>

#include "PDE_SWESphere2D/BenchmarksCombined.hpp"

#include <sweet/Error/Fatal.hpp>


/*
 * This allows running REXI including Coriolis-related terms but just by setting f to 0
 */


class ProgramVisSphericalHarmonics
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

		sweet::Data::Sphere2D::DataSpectral sphere2DDataForVisualization;


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
				sweet::Data::Sphere2D::Shack *i_shackSphere2DDataOps
		)
		{
			/*
			 * Setup Sphere2D Data Config & Operators
			 */
			sphere2DDataConfig.setupAuto(i_shackSphere2DDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphere2DDataConfig);

			ops.setup(&sphere2DDataConfig, i_shackSphere2DDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(ops);

			sphere2DDataForVisualization.setup(sphere2DDataConfig);

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
			sphere2DDataForVisualization.clear();

			ops.clear();
			sphere2DDataConfig.clear();
		}
	};

	// Simulation data
	DataConfigOps dataConfigOps;

	// Diagnostics measures
	int last_timestep_nr_update_diagnostics = -1;


#if SWEET_GUI && 0
	// Data to visualize is stored to this variable
	sweet::Data::Cart2D::DataGrid vis_cart2d_data;

	// Which primitive to use for rendering
	int vis_render_type_of_primitive_id = 1;

	// Which primitive to use for rendering
	int vis_data_id = 0;
#endif


	// was the output of the time step already done for this simulation state?
	double timestep_last_output_simtime;

	int mode_m = 0;
	int mode_n = 0;

	bool vis_reset = false;


	/*
	 * Shack directory and shacks to work with
	 */
	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict;
	sweet::Data::Sphere2D::Shack *shackSphere2DDataOps;
	sweet::IO::Shack *shackIOData;



public:
	ProgramVisSphericalHarmonics(
			int i_argc,
			char *const * const i_argv
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackSphere2DDataOps(nullptr),
		shackIOData(nullptr)
	{
	}

	~ProgramVisSphericalHarmonics()
	{
		clear();
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
		shackSphere2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Sphere2D::Shack>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);
		shackIOData = shackProgArgDict.getAutoRegistration<sweet::IO::Shack>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

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
		shackIOData = nullptr;

		shackProgArgDict.clear();
	}

	bool setup_2_processArguments()
	{
		/*
		 * SHACK: Process arguments
		 */
		shackProgArgDict.processProgramArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		return true;
	}

	void clear_2_processArguments()
	{
		shackProgArgDict.clear();
	}

	bool setup_3_data()
	{
		/*
		 * Setup the data fields
		 */
		dataConfigOps.setup(shackSphere2DDataOps);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataConfigOps);

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


	void setup_mode()
	{
		if (mode_n < 0)
			mode_n = 0;

		if (mode_m < 0)
			mode_m = 0;

		if (mode_n > dataConfigOps.sphere2DDataConfig.spectral_modes_n_max)
			mode_n = dataConfigOps.sphere2DDataConfig.spectral_modes_n_max;

		if (mode_m > dataConfigOps.sphere2DDataConfig.spectral_modes_m_max)
			mode_m = dataConfigOps.sphere2DDataConfig.spectral_modes_m_max;

		if (mode_m > mode_n)
			mode_m = mode_n;

		dataConfigOps.sphere2DDataForVisualization.spectral_setZero();
		std::complex<double> val = 1;
		dataConfigOps.sphere2DDataForVisualization.spectral_set(mode_n, mode_m, val);
	}

public:
	bool should_quit()
	{
		return false;
	}


	bool runTimestep()
	{
		return true;
	}



#if SWEET_GUI
	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(
			int i_num_iterations
	)
	{
		runTimestep();
	}


	int max_vis_types = 9;


	void vis_getDataArray(
			const sweet::Data::Cart2D::DataGrid **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive_id,
			void **o_bogus_data,
			double *o_vis_min,
			double *o_vis_max,
			bool *o_vis_reset
	)
	{
		*o_vis_reset = vis_reset;
		vis_reset = false;

		// request rendering of sphere2D or cart2d
		*o_render_primitive_id = dataConfigOps.vis_render_type_of_primitive_id;
		*o_bogus_data = &dataConfigOps.sphere2DDataConfig;


		//int id = dataConfigOps.vis_data_id % max_vis_types;

		dataConfigOps.vis_cart2d_data = sweet::Data::Sphere2D::Convert::DataSpectral_2_Cart2D_DataGrid::convert(dataConfigOps.sphere2DDataForVisualization, dataConfigOps.cart2DDataConfig);

		double vis_min = dataConfigOps.vis_cart2d_data.grid_reduce_min();
		double vis_max = dataConfigOps.vis_cart2d_data.grid_reduce_max();

		vis_max = std::max(std::abs(vis_max), std::abs(vis_min));
		vis_min = -vis_max;

		*o_vis_min = vis_min;
		*o_vis_max = vis_max;


		*o_dataArray = &dataConfigOps.vis_cart2d_data;
		*o_aspect_ratio = 0.5;
	}



	/**
	 * return status string for window title
	 */
	const std::string vis_getStatusString(bool &o_replace_commas_with_newline)
	{
		o_replace_commas_with_newline = true;
		std::ostringstream description;

		description << "Visualizing modes (N=" << mode_n << " | M=" << mode_m << ")";
		description << ",";
		description << "  Min: " << dataConfigOps.vis_cart2d_data.grid_reduce_max();
		description << ",";
		description << "  Max: " << dataConfigOps.vis_cart2d_data.grid_reduce_max();

		return description.str();
	}



	void vis_pause()
	{
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
			dataConfigOps.vis_render_type_of_primitive_id = (dataConfigOps.vis_render_type_of_primitive_id + 1) % 2;
			break;

		case 'm':
			mode_m++;
			setup_mode();
			break;

		case 'M':
			mode_m--;
			setup_mode();
			break;

		case 'n':
			mode_n++;
			setup_mode();
			break;

		case 'N':
			mode_n--;
			setup_mode();
			break;

		case 't':
			shackSphere2DDataOps->space_res_spectral[0] *= 2;
			shackSphere2DDataOps->space_res_spectral[1] *= 2;
			reset();
			break;

		case 'T':
			shackSphere2DDataOps->space_res_spectral[0] /= 2;
			shackSphere2DDataOps->space_res_spectral[1] /= 2;
			reset();
			break;
		}
	}
#endif
};



int main(int i_argc, char *i_argv[])
{
	ProgramVisSphericalHarmonics visSphericalHarmonics(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(visSphericalHarmonics);

	visSphericalHarmonics.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(visSphericalHarmonics);

#if SWEET_GUI
	if (visSphericalHarmonics.shackIOData->guiEnabled)
	{
		sweet::GUI::VisSweet visSweet(visSphericalHarmonics);
	}
	else
#endif
	{
		while (!visSphericalHarmonics.should_quit())
			visSphericalHarmonics.runTimestep();
	}

	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(visSphericalHarmonics);


	std::cout << "FIN" << std::endl;
	return 0;
}
