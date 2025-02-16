/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_AdvectionSphere2D/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_AdvectionSphere2D/time/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_AdvectionSphere2D/benchmarks/
 *
 * MULE_SCONS_OPTIONS: --sphere2d-spectral-space=enable
 */


#include <sweet/Data/Sphere2D/Convert/DataGrid_2_Cart2D_DataGrid.hpp>
#include <sweet/Data/Sphere2D/Convert/DataSpectral_2_Cart2D_DataGrid.hpp>
#include <programs/PDE_AdvectionSphere2D/PDEAdvectionSphere2DBenchmarksCombined.hpp>
#include <programs/PDE_AdvectionSphere2D/PDEAdvectionSphere2DTimeSteppers.hpp>
#include <programs/PDE_AdvectionSphere2D/ProgramPDEAdvectionSphere2D.hpp>
#include <sweet/Data/Sphere2D/DataSpectral.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/SemiLagrangian/Shack.hpp>
#include <sweet/SemiLagrangian/Sphere2D.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include "../programs/PDE_AdvectionSphere2D/time/ShackPDEAdvectionSphere2DTimeDiscretization.hpp"




int main(int i_argc, char *i_argv[])
{
	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict(i_argc, i_argv);
	shackProgArgDict.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	sweet::Data::Sphere2D::Shack *shackSphere2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Sphere2D::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	sweet::TimeTree::Shack *shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::TimeTree::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	sweet::SemiLagrangian::Shack *shackTimesteppingSemiLagrangianSphere2DData = shackProgArgDict.getAutoRegistration<sweet::SemiLagrangian::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	ShackPDEAdvectionSphere2DTimeDiscretization *shackTimeDisc = shackProgArgDict.getAutoRegistration<ShackPDEAdvectionSphere2DTimeDiscretization>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.printShackData();

	int initial_spectral_modes = shackSphere2DDataOps->space_res_spectral[0];

	if (shackTimestepControl->currentTimestepSize < 0)
		SWEETErrorFatal("Timestep size not set");

	int max_modes = 256;

	if (shackTimeDisc->timestepping_order == 1)
		max_modes = 512;
	else if (shackTimeDisc->timestepping_order == 2)
		max_modes = 256;

	double prev_max_error = -1;
	for (int i = initial_spectral_modes; i <= max_modes; i *= 2)
	{
		shackTimestepControl->currentTimestepSize *= 0.5;
		//shackTimestepControl->setup_timestepSize = shackTimestepControl->current_timestepSize;

		if (shackTimeDisc->timestepping_method == "na_sl")
		{
			shackSphere2DDataOps->space_res_spectral[0] = i;
			shackSphere2DDataOps->space_res_spectral[1] = i;

			shackSphere2DDataOps->space_res_physical[0] = 2*i;
			shackSphere2DDataOps->space_res_physical[1] = i;
		}
		else
		{
			shackSphere2DDataOps->space_res_spectral[0] = initial_spectral_modes;
			shackSphere2DDataOps->space_res_spectral[1] = initial_spectral_modes;

			shackSphere2DDataOps->space_res_physical[0] = 0;
			shackSphere2DDataOps->space_res_physical[1] = 0;
		}


		if (1)
		{
			sweet::Data::Sphere2D::Config sphere2DDataConfigInstance;
			sphere2DDataConfigInstance.setupAuto(
					shackSphere2DDataOps
			);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(sphere2DDataConfigInstance);

			sweet::Data::Sphere2D::Config *sphere2DDataConfig = &sphere2DDataConfigInstance;

			std::cout << "Checking for right velocity" << std::endl;
			sweet::SemiLagrangian::Sphere2D sl;
			sl.setup(sphere2DDataConfig, shackTimesteppingSemiLagrangianSphere2DData, shackTimeDisc->timestepping_order);


			/*
			 * Convert to Cartesian velocity space
			 */
			sweet::Data::Vector::Vector<double> V_lon_D(sphere2DDataConfig->grid_number_elements);
			sweet::Data::Vector::Vector<double> V_lat_D(sphere2DDataConfig->grid_number_elements);
			V_lon_D.set_all(1.0);
			V_lat_D.set_all(2.0);

			sweet::Data::Vector::Vector<double> V_x_D, V_y_D, V_z_D;
			sweet::LibMath::VectorMath::velocity_latlon_2_cartesian__array(
					sl.pos_lon_A,
					sl.pos_lat_A,
					V_lon_D,
					V_lat_D,
					V_x_D,
					V_y_D,
					V_z_D
				);

			sweet::Data::Vector::Vector<double> V_lon_tmp, V_lat_tmp;
			sweet::LibMath::VectorMath::velocity_cartesian_2_latlon__array(
					sl.pos_lon_A,
					sl.pos_lat_A,
					V_x_D,
					V_y_D,
					V_z_D,
					V_lon_tmp, V_lat_tmp
			);

			double err_lon = (V_lon_D - V_lon_tmp).reduce_maxAbs();
			double err_lat = (V_lat_D - V_lat_tmp).reduce_maxAbs();

			if (err_lon > 1e-10)
			{
				std::cerr << "Error: " << err_lon << std::endl;
				SWEETErrorFatal("Error lon too high!");
			}

			if (err_lat > 1e-10)
			{
				std::cerr << "Error: " << err_lat << std::endl;
				SWEETErrorFatal("Error lat too high!");
			}
		}

		std::cout << "Testing with Sphere2DDataOps:" << std::endl;
		shackSphere2DDataOps->printShack();
		std::cout << std::endl;

		std::cout << "Testing with dt=" << shackTimestepControl->currentTimestepSize << std::endl;

		ProgramPDEAdvectionSphere2D simulation(i_argc, i_argv);

		simulation.setup_1_shackRegistration();
		simulation.setup_2_processArguments();

		// physical resolution
		simulation.shackSphere2DDataOps->space_res_physical[0] = shackSphere2DDataOps->space_res_physical[0];
		simulation.shackSphere2DDataOps->space_res_physical[1] = shackSphere2DDataOps->space_res_physical[1];

		// spectral resolution
		simulation.shackSphere2DDataOps->space_res_spectral[0] = shackSphere2DDataOps->space_res_spectral[0];
		simulation.shackSphere2DDataOps->space_res_spectral[1] = shackSphere2DDataOps->space_res_spectral[1];

		// time step size
		simulation.shackTimestepControl->currentTimestepSize = shackTimestepControl->currentTimestepSize;

		// Setup the data and operators
		simulation.setup_3_dataOpsEtc();

		{
			while (!simulation.should_quit())
				simulation.runTimestep();

			double error_lmax = simulation.getErrorLMaxOnH();
			double error_rms = simulation.getErrorRMSOnH();

			std::cout << "Error compared to initial condition" << std::endl;
			std::cout << "Lmax error: " << error_lmax << std::endl;
			std::cout << "RMS error: " << error_lmax << std::endl;

			if (prev_max_error >= 0)
			{
				//double conv = (prev_max_error - simulation.max_error) / simulation.max_error;
				double conv = prev_max_error / error_lmax;
				std::cout << "Convergence: " << conv << std::endl;

				if (conv*1.1 < std::pow(2.0, (double)shackTimeDisc->timestepping_order))
					SWEETErrorFatal("Convergence not given!");
			}

			if (error_lmax  > 1e10)
				SWEETErrorFatal("Lmax error exceeded threshold!");

			prev_max_error = error_lmax;

			std::cout << "*********************************************" << std::endl;
		}
	}

	return 0;
}
