/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_AdvectionSphere2D/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_AdvectionSphere2D/time/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_AdvectionSphere2D/benchmarks/
 *
 * MULE_SCONS_OPTIONS: --sphere2d-spectral-space=enable
 */

#include <programs/PDE_AdvectionSphere2D/ProgramPDEAdvectionSphere2D.hpp>

int main(int i_argc, char *i_argv[])
{
	ProgramPDEAdvectionSphere2D progPDEAdvSphere2D(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(progPDEAdvSphere2D);

	progPDEAdvSphere2D.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(progPDEAdvSphere2D);

	// Simply test whether the clear and setup works properly
	progPDEAdvSphere2D.clear();
	progPDEAdvSphere2D.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(progPDEAdvSphere2D);


	int max_modes = 256;

	double dt = progPDEAdvSphere2D.shackTimestepControl->currentTimestepSize;
	int initial_spectral_modes = progPDEAdvSphere2D.shackSphere2DDataOps->space_res_spectral[0];

	int loop_counter = 0;

	double prev_max_error = -1;
	for (int i = initial_spectral_modes; i <= max_modes; i *= 2)
	{
		// only clear and setup the spatial parts
		progPDEAdvSphere2D.clear();

		// setup stage 1 + 2
		progPDEAdvSphere2D.setup_1_shackRegistration();
		progPDEAdvSphere2D.setup_2_processArguments();

		/*
		 * Overwrite parameters
		 */
		progPDEAdvSphere2D.shackTimestepControl->currentTimestepSize = dt/std::pow(2.0, loop_counter);

		/*
		 * We need higher resolution for
		 * - SL methods or if using
		 * - finite differences in space
		 */
		int expected_order;

		if (
				progPDEAdvSphere2D.shackTimeDisc->timestepping_method == "na_sl"
		)
		{
			progPDEAdvSphere2D.shackSphere2DDataOps->space_res_spectral[0] = i;
			progPDEAdvSphere2D.shackSphere2DDataOps->space_res_spectral[1] = i;

			progPDEAdvSphere2D.shackSphere2DDataOps->space_res_physical[0] = 0;
			progPDEAdvSphere2D.shackSphere2DDataOps->space_res_physical[1] = 0;

			// order is limited by the spatial interpolation order
			expected_order = std::min(progPDEAdvSphere2D.shackTimeDisc->timestepping_order, 2);
		}
		else
		{
			progPDEAdvSphere2D.shackSphere2DDataOps->space_res_spectral[0] = initial_spectral_modes;
			progPDEAdvSphere2D.shackSphere2DDataOps->space_res_spectral[1] = initial_spectral_modes;

			progPDEAdvSphere2D.shackSphere2DDataOps->space_res_physical[0] = 0;
			progPDEAdvSphere2D.shackSphere2DDataOps->space_res_physical[1] = 0;

			// order is directly the time discretization order
			expected_order = progPDEAdvSphere2D.shackTimeDisc->timestepping_order;
		}

		// setup 3rd part
		progPDEAdvSphere2D.setup_3_dataOpsEtc();

		std::cout << "Testing with " << progPDEAdvSphere2D.dataConfigOps.sphere2DDataConfig.getUniqueIDString() << std::endl;
		std::cout << "Testing with dt=" << progPDEAdvSphere2D.shackTimestepControl->currentTimestepSize << std::endl;

		{
			while (!progPDEAdvSphere2D.should_quit())
				progPDEAdvSphere2D.runTimestep();

			double max_error = progPDEAdvSphere2D.getErrorLMaxOnH();

			std::cout << "Lmax error compared to initial condition: " << max_error << std::endl;

			if (prev_max_error >= 0)
			{
				//double conv = (prev_max_error - simulation.max_error) / simulation.max_error;
				double conv = prev_max_error / max_error;
				std::cout << "prev_max_error: " << prev_max_error << std::endl;
				std::cout << "max_error: " << max_error << std::endl;
				std::cout << "Convergence: " << conv << std::endl;
				std::cout << "expected_order: " << expected_order << std::endl;

				// TODO: We are much more tolerant here
				if (conv < std::pow(2.0, (double)expected_order)*0.7)
				{
					std::cerr << "Convergence too low!" << std::endl;
					exit(1);
				}

				// TODO: We need to be very tolerant with a too high order (2.0 instead of 1.3)
				if (conv > std::pow(2.0, (double)(expected_order))*2.0)
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
