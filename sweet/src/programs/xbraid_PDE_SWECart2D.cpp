/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWECart2D/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWECart2D/DataContainer/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWECart2D/TimeTree/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWECart2D/TimeTree/TimeStepper
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWECart2D/Benchmarks/
 *
 *
 * MULE_SCONS_OPTIONS: --cart2d-spectral-space=enable
 * MULE_SCONS_OPTIONS: --cart2d-spectral-dealiasing=enable
 * MULE_SCONS_OPTIONS: --sweet-mpi=enable
 * MULE_SCONS_OPTIONS: --xbraid=mpi
 * MULE_SCONS_OPTIONS: --xbraid-cart2d=enable
 * MULE_SCONS_OPTIONS: --xbraid-cart2d-swe=enable
 */

#include <programs/PDE_SWECart2D/XBraid/Program.hpp>

int main(int i_argc, char *i_argv[])
{

#if SWEET_MPI
	int mpi_rank;
	MPI_Init(&i_argc, &i_argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

	PDE_SWECart2D::XBraid::Program simulation(
						i_argc, i_argv
#if SWEET_MPI
						,
						MPI_COMM_WORLD,
						mpi_rank
#endif
			);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	simulation.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	{
		simulation.shackTimestepControl->validateMaxSimulationTimeOrTimestepNr();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*(simulation.shackTimestepControl));

		simulation.runXBraid();

	}

	///simulation.printSimulationErrors();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	std::cout << "FIN" << std::endl;

#if SWEET_MPI
	MPI_Finalize();
#endif

	return 0;
}
