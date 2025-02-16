/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 * 
 * MULE_SCONS_OPTIONS: --cart2d-spectral-space=enable
 * MULE_SCONS_OPTIONS: --gui=enable
 */

#include <sweet/Tools/DefaultPrecompilerValues.hpp>
#include <ostream>
#include <sstream>

// This is just for the editor to show code as used within precompiler #if ... directives

#if !SWEET_USE_CART2D_SPECTRAL_SPACE
	#error "Spectral space required"
#endif

// Error handling
#include <sweet/Error/Base.hpp>

// Parse program arguments
#include <sweet/Tools/ProgramArguments.hpp>

// Include everything we need for simulations on the cart2d
#include <sweet/Data/Cart2D/Cart2D.hpp>

// Our shack directory to store different objects and get them back later on
#include <sweet/Shacks/Dictionary.hpp>

#include <sweet/Tools/Stopwatch.hpp>

#include <sweet/GUI/VisSweet.hpp>



class ProgramCart2DSpectralVisualization
#if SWEET_GUI
		:	public sweet::GUI::SimulationGUICallbacks
#endif
{
public:
	sweet::Error::Base error;
	sweet::Tools::ProgramArguments programArguments;

	class Data
	{
	public:
		sweet::Error::Base error;

		sweet::Data::Cart2D::Config cart2DDataConfig;
		sweet::Data::Cart2D::Operators ops;

		sweet::Data::Cart2D::DataSpectral tmp;
		sweet::Data::Cart2D::DataGrid tmp_phys;

		bool setup(sweet::Data::Cart2D::Shack *i_shackCart2DDataOps)
		{
			/*
			 * Setup Cart2D Data Config & Operators
			 */
			if (!cart2DDataConfig.setupAuto(*i_shackCart2DDataOps))
				return error.forwardWithPositiveReturn(cart2DDataConfig.error);

			if (!ops.setup(
					cart2DDataConfig,
					*i_shackCart2DDataOps
				))
				return error.forwardWithPositiveReturn(ops.error);

			tmp.setup(cart2DDataConfig);
			tmp_phys.setup(cart2DDataConfig);

			return true;
		}

		void clear()
		{
			tmp.clear();
			tmp_phys.clear();

			ops.clear();
			cart2DDataConfig.clear();
		}
	};

	// Data
	Data data;
	
	sweet::Shacks::Dictionary shackDict;
	sweet::Data::Cart2D::Shack *shackCart2DDataOps;


	// Data to visualize is stored to this variable
	sweet::Data::Cart2D::DataGrid vis_cart2d_data;

	// Which primitive to use for rendering
	int vis_dataId = 0;

	std::string vis_description;


	sweet::Tools::Stopwatch stopwatch;

	/*
	 * Setup main
	 */
	int prog_argc;
	char *const * const prog_argv;

public:
	ProgramCart2DSpectralVisualization(
			int i_argc,
			char *const * const i_argv
	)	:
		prog_argc(i_argc),
		prog_argv(i_argv)
	{
	}


	/*
	 * Clear all data
	 */
	void clear()
	{
		programArguments.clear();
		shackDict.clear();

		data.clear();

#if SWEET_GUI
		vis_cart2d_data.clear();
#endif
	}

	/*
	 * Chekc if help should be printed a do so
	 */
	void checkAndPrintHelp()
	{
		/*
		 * First, check for --help or -h
		 */
		if (programArguments.argumentWithKeyExists("-h") || programArguments.argumentWithKeyExists("--help"))
		{
			std::cout << "Printing help:" << std::endl;
			shackDict.printProgramArguments();
		}
	}


	bool setup()
	{
		/*
		 * Parse program arguments
		 */
		programArguments.setup(prog_argc, prog_argv);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(programArguments);

		/*
		 * SHACK: Register classes which we require
		 */
		shackCart2DDataOps = shackDict.getAutoRegistration<sweet::Data::Cart2D::Shack>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackDict);

		/*
		 * First, check for --help or -h
		 */
		if (programArguments.argumentWithKeyExists("-h") || programArguments.argumentWithKeyExists("--help"))
		{
			std::cout << "Printing help:" << std::endl;
			shackDict.printProgramArguments();
			return false;
		}

		/*
		 * SHACK: Process arguments
		 */
		shackDict.processProgramArguments(programArguments);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackDict);

		/*
		 * Setup Cart2D Data Config & Operators
		 */
		data.setup(shackCart2DDataOps);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(data);


#if SWEET_GUI
		vis_cart2d_data.setup(data.cart2DDataConfig);
#endif

		std::cout << "Printing shack information:" << std::endl;
		shackDict.printShackData();

		/*
		 * Finish registration & getting class interfaces so that nobody can do some
		 * strange things with this anymore
		 */
		shackDict.closeRegistration();
		shackDict.closeGet();

		/*
		 * Now we should check that all program arguments have really been parsed
		 */
		programArguments.checkAllArgumentsProcessed();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(programArguments);

		return true;
	}

	bool reset()
	{
		clear();
		setup();

		return !error.exists();
	}


	bool runTimestep()
	{
		return true;
	}


public:
	void timestep_output(
			std::ostream &o_ostream = std::cout
	)
	{
	}



public:
	bool should_quit()
	{
		return false;
	}


	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(
		int i_num_iterations
	)
	{
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
		data.tmp.spectral_setZero();
		int spec_array[][4] =
		{
				// j, i, re, im
				{0, 0, 1, 0},
				{0, 0, 0, 1},
				{0, 0, 1, 1},

				{1, 0, 1, 0},
				{1, 0, 0, 1},
				{1, 0, 1, 1},

				{2, 0, 1, 0},
				{2, 0, 0, 1},
				{2, 0, 1, 1},

				{0, 1, 1, 0},
				{0, 1, 0, 1},
				{0, 1, 1, 1},

				{1, 1, 1, 0},
				{1, 1, 0, 1},
				{1, 1, 1, 1},

				{(int)shackCart2DDataOps->space_res_physical[1]-1, 0, 1, 0},
				{(int)shackCart2DDataOps->space_res_physical[1]-1, 0, 0, 1},
				{(int)shackCart2DDataOps->space_res_physical[1]-1, 0, 1, 1},

				{(int)shackCart2DDataOps->space_res_physical[1]-2, 0, 1, 0},
				{(int)shackCart2DDataOps->space_res_physical[1]-2, 0, 0, 1},
				{(int)shackCart2DDataOps->space_res_physical[1]-2, 0, 1, 1},

				{0, (int)shackCart2DDataOps->space_res_physical[0]/2, 1, 0},
				{0, (int)shackCart2DDataOps->space_res_physical[0]/2, 0, 1},
				{0, (int)shackCart2DDataOps->space_res_physical[0]/2, 1, 1},

				{0, (int)shackCart2DDataOps->space_res_physical[0]/2-1, 1, 0},
				{0, (int)shackCart2DDataOps->space_res_physical[0]/2-1, 0, 1},
				{0, (int)shackCart2DDataOps->space_res_physical[0]/2-1, 1, 1},

				{0, (int)shackCart2DDataOps->space_res_physical[0]/2-2, 1, 0},
				{0, (int)shackCart2DDataOps->space_res_physical[0]/2-2, 0, 1},
				{0, (int)shackCart2DDataOps->space_res_physical[0]/2-2, 1, 1},
		};

		int id = vis_dataId % (sizeof(spec_array)/sizeof(spec_array[0]));

		std::ostringstream ss;
		ss << "spec_coord (j, i) = (" << spec_array[id][0] << ", " << spec_array[id][1] << ")";
		ss << ", ";
		ss << "value = (" << spec_array[id][2] << ", " << spec_array[id][3] << "i)";
		vis_description = ss.str();

		data.tmp.spectral_set(
				spec_array[id][0],
				spec_array[id][1],
				{
						(double)spec_array[id][2],
						(double)spec_array[id][3]
				}
			);

		data.tmp_phys = data.tmp.toGrid();

		*o_dataArray = &data.tmp_phys;

		*o_aspect_ratio = shackCart2DDataOps->cart2d_domain_size[1] / shackCart2DDataOps->cart2d_domain_size[0];
	}


	/**
	 * return status string for window title
	 */
	const std::string vis_getStatusString(bool &o_replace_commas_with_newline)
	{
		o_replace_commas_with_newline = false;
		return vis_description;
	}



	void vis_pause()
	{
	}



	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			vis_dataId++;
			break;

		case 'V':
			vis_dataId--;
			break;
		}
	}


};


int main(int i_argc, char *i_argv[])
{
	ProgramCart2DSpectralVisualization simulation(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	simulation.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);


	sweet::GUI::VisSweet visSweet(simulation);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	std::cout << "FIN" << std::endl;
	return 0;
}
