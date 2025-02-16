/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/ODE_Generic/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/ODE_Generic/TimeTree
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/ODE_Generic/DE_Dahlquist/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/ODE_Generic/DE_Dahlquist/TimeTree
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/ODE_Generic/DE_Dahlquist/TimeTree/DE_Terms
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/ODE_Generic/DE_Dahlquist/Benchmarks
 */

#if SWEET_GUI
#	error "GUI not supported for this program"
#endif

#include <sweet/Tools/DefaultPrecompilerValues.hpp>

#if SWEET_MPI
#	include <mpi.h>
#endif

#include "ODE_Generic/Program.hpp"


#if SWEET_MPI
int mpi_comm_rank;
int mpi_comm_size;
#endif

bool isMPIRoot()
{
#if SWEET_MPI
	return mpi_comm_rank == 0;
#else
	return true;
#endif
}


int main_mpi(int i_argc, char *i_argv[])
{
#if SWEET_MPI
	#if SWEET_THREADING_SPACE
		int provided;
		MPI_Init_thread(&i_argc, &i_argv, MPI_THREAD_MULTIPLE, &provided);

		if (provided != MPI_THREAD_MULTIPLE)
			SWEETErrorFatal("MPI_THREAD_MULTIPLE not available! Try to get an MPI version with multi-threading support or compile without OMP/TBB support. Good bye...");
	#else
		MPI_Init(&i_argc, &i_argv);
	#endif

	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_comm_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_comm_size);
#endif

	ODE_Generic::Program simulation(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	simulation.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);


	simulation.shackTimeTree->validateMaxSimulationTimeOrTimestepNr();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*(simulation.shackTimeTree));

	simulation.timestepHandleOutput();

	sweet::Tools::StopwatchBox::getInstance().main_timestepping.start();

	while (!simulation.should_quit())
	{
		simulation.runTimestep();

		if (isMPIRoot())
			simulation.timestepHandleOutput();
	}

	if (isMPIRoot())
		std::cout << "TIMESTEPPING FINISHED" << std::endl;

	sweet::Tools::StopwatchBox::getInstance().main_timestepping.stop();


	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	if (isMPIRoot())
		std::cout << "FIN" << std::endl;

#if SWEET_MPI
	MPI_Finalize();
#endif

	return 0;
}


int main(int i_argc, char *i_argv[])
{
	int retval = main_mpi(i_argc, i_argv);

	return retval;
}

