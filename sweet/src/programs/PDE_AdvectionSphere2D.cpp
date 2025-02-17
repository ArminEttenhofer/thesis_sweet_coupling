/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_AdvectionSphere2D/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_AdvectionSphere2D/time
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_AdvectionSphere2D/benchmarks
 *
 * MULE_SCONS_OPTIONS: --sphere2d-spectral-space=enable
 */

#ifndef SWEET_GUI
	#define SWEET_GUI 1
#endif

#include <iostream>

#include <programs/PDE_AdvectionSphere2D/ProgramPDEAdvectionSphere2D.hpp>


int main(int i_argc, char *i_argv[])
{
	ProgramPDEAdvectionSphere2D simulation(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	simulation.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

#if SWEET_GUI
	if (simulation.shackIOData->guiEnabled)
	{
		sweet::GUI::VisSweet visSweet(simulation);
	}
	else
#endif
	{
		simulation.shackTimestepControl->validateMaxSimulationTimeOrTimestepNr();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*(simulation.shackTimestepControl));

		while (!simulation.should_quit())
			simulation.runTimestep();
	}

	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	if (simulation.shackPDEAdvectionSphere2D->compute_errors)
	{
		simulation.printSimulationErrors();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);
	}

	std::cout << "FIN" << std::endl;
	return 0;
}
