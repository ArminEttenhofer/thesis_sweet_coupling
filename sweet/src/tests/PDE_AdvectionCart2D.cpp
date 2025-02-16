/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_AdvectionCart2D/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_AdvectionCart2D/time/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_AdvectionCart2D/benchmarks/
 *
 * MULE_SCONS_OPTIONS: --cart2d-spectral-space=enable
 */

#include <programs/PDE_AdvectionCart2D/ProgramPDEAdvectionCart2D.hpp>

int main(int i_argc, char *i_argv[])
{
	ProgramPDEAdvectionCart2D progPDEAdvCart2D(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(progPDEAdvCart2D);

	progPDEAdvCart2D.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(progPDEAdvCart2D);

	// Simply test whether the clear and setup works properly
	progPDEAdvCart2D.clear();
	progPDEAdvCart2D.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(progPDEAdvCart2D);


	int max_modes = 256;

	double dt = progPDEAdvCart2D.shackTimestepControl->currentTimestepSize;
	int initial_spectral_modes = progPDEAdvCart2D.shackCart2DDataOps->space_res_spectral[0];

	int loop_counter = 0;

	double prev_max_error = -1;
	for (int i = initial_spectral_modes; i <= max_modes; i *= 2)
	{
		// only clear and setup the spatial parts
		progPDEAdvCart2D.clear();

		// setup stage 1 + 2
		progPDEAdvCart2D.setup_1_shackRegistration();
		progPDEAdvCart2D.setup_2_processArguments();

		/*
		 * Overwrite parameters
		 */
		progPDEAdvCart2D.shackTimestepControl->currentTimestepSize = dt/std::pow(2.0, loop_counter);

		/*
		 * We need higher resolution for
		 * - SL methods or if using
		 * - finite differences in space
		 */
		int expected_order;

		if (
				progPDEAdvCart2D.shackTimeDisc->timestepping_method == "na_sl" ||
				progPDEAdvCart2D.shackCart2DDataOps->space_use_spectral_basis_diffs == false
		)
		{
			progPDEAdvCart2D.shackCart2DDataOps->space_res_spectral[0] = i;
			progPDEAdvCart2D.shackCart2DDataOps->space_res_spectral[1] = i;

			progPDEAdvCart2D.shackCart2DDataOps->space_res_physical[0] = 0;
			progPDEAdvCart2D.shackCart2DDataOps->space_res_physical[1] = 0;

			// order is limited by the spatial interpolation order
			expected_order = std::min(progPDEAdvCart2D.shackTimeDisc->timestepping_order, 2);
		}
		else
		{
			progPDEAdvCart2D.shackCart2DDataOps->space_res_spectral[0] = initial_spectral_modes;
			progPDEAdvCart2D.shackCart2DDataOps->space_res_spectral[1] = initial_spectral_modes;

			progPDEAdvCart2D.shackCart2DDataOps->space_res_physical[0] = 0;
			progPDEAdvCart2D.shackCart2DDataOps->space_res_physical[1] = 0;

			// order is directly the time discretization order
			expected_order = progPDEAdvCart2D.shackTimeDisc->timestepping_order;
		}

		// setup 3rd part
		progPDEAdvCart2D.setup_3_dataOpsEtc();

		std::cout << "Testing with " << progPDEAdvCart2D.dataConfigOps.cart2DDataConfig.getUniqueIDString() << std::endl;
		std::cout << "Testing with dt=" << progPDEAdvCart2D.shackTimestepControl->currentTimestepSize << std::endl;

		{
			while (!progPDEAdvCart2D.should_quit())
				progPDEAdvCart2D.runTimestep();

			double max_error = progPDEAdvCart2D.getErrorLMaxOnH();

			std::cout << "Lmax error compared to initial condition: " << max_error << std::endl;

			if (prev_max_error >= 0)
			{
				//double conv = (prev_max_error - simulation.max_error) / simulation.max_error;
				double conv = prev_max_error / max_error;
				std::cout << "prev_max_error: " << prev_max_error << std::endl;
				std::cout << "max_error: " << max_error << std::endl;
				std::cout << "Convergence: " << conv << std::endl;
				std::cout << "expected_order: " << expected_order << std::endl;

				if (conv < std::pow(2.0, (double)expected_order)*0.9)
				{
					std::cerr << "Convergence too low!" << std::endl;
					exit(1);
				}

				if (conv > std::pow(2.0, (double)(expected_order))*1.1)
				{
					std::cerr << "Convergence too high, stopping here!" << std::endl;
					exit(1);
				}
			}
			prev_max_error = max_error;
			std::cout << std::endl;
		}

		loop_counter++;
	}

	return 0;
}
