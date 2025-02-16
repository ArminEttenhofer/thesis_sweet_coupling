/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_AdvectionCart2D/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_AdvectionCart2D/time/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_AdvectionCart2D/benchmarks/
 *
 * MULE_SCONS_OPTIONS: --cart2d-spectral-space=enable
 */

#ifndef SWEET_GUI
	#define SWEET_GUI 1
#endif

#include <programs/PDE_AdvectionCart2D/ProgramPDEAdvectionCart2D.hpp>


int main(int i_argc, char *i_argv[])
{
	ProgramPDEAdvectionCart2D simulation(i_argc, i_argv);
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

	simulation.printSimulationErrors();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	std::cout << "FIN" << std::endl;
	return 0;
}
