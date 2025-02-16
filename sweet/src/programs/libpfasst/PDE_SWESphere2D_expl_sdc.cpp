/*
 * Author: Valentina Schüller & Francois Hamon & Martin Schreiber <SchreiberX@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/libpfasst/PDE_SWESphere2D_expl_sdc/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWESphere2D/TimeOld/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWESphere2D/Benchmarks/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWESphere2D/
 *
 * MULE_SCONS_OPTIONS: --sphere2d-spectral-space=enable
 * MULE_SCONS_OPTIONS: --fortran-source=enable
 * MULE_SCONS_OPTIONS: --sweet-mpi=enable
 * MULE_SCONS_OPTIONS: --libpfasst=enable
 */

#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <iostream>
#include <mpi.h>
#include <sweet/Error/Fatal.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>

#include "../PDE_SWESphere2D/Shack.hpp"
#include "../PDE_SWESphere2D/Benchmarks/Shack.hpp"
#include "ShackLibPFASST.hpp"

#include "interface/LevelSingleton.hpp"
#include "PDE_SWESphere2D_expl_sdc/Sphere2DDataCtxSDC.hpp"
#include "../PDE_SWESphere2D/BenchmarksCombined.hpp"
#include "../PDE_SWESphere2D/Diagnostics.hpp"
#include "../PDE_SWESphere2D/TimeOld/ShackTimeDiscretization.hpp"

#define WITH_MPI

extern "C"
{
	/* Driver function for pfasst control */
	void fmain(Sphere2DDataCtxSDC *pd_ctx,
			   const int *niters,
			   const int *nsweeps_coarse,
			   const int *nnodes,
			   const char *qtype_name,
			   const int *qtype_name_len,
			   const int *use_rk_stepper, // 1 means true, 0 means false
			   const int *nfields,
			   const int *nvars_per_field,
			   double *t_max,
			   double *dt);
}

/**
 * Main function launching LibPFASST
 */

int main(int i_argc, char *i_argv[])
{
	MPI_Init(&i_argc, &i_argv);

	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict(i_argc, i_argv);
	shackProgArgDict.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	sweet::Data::Sphere2D::Shack *shackSphere2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Sphere2D::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	//PDE_SWESphere2D::Shack *shackPDESWESphere2D =
			shackProgArgDict.getAutoRegistration<PDE_SWESphere2D::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	ShackLibPFASST *shackLibPFASST = shackProgArgDict.getAutoRegistration<ShackLibPFASST>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	// sweet::TimeTree::ShackTimestepControl *shackTimestepControl =
	shackProgArgDict.getAutoRegistration<sweet::TimeTree::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	//PDE_SWESphere2D::Benchmarks::Shack *shackBenchmarks =
		shackProgArgDict.getAutoRegistration<PDE_SWESphere2D::Benchmarks::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	// sweet::IO::ShackIOData *shackIOData =
	shackProgArgDict.getAutoRegistration<sweet::IO::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	// ShackTimeDiscretization *shackTimeDisc =
	shackProgArgDict.getAutoRegistration<ShackTimeDiscretization>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.checkAllArgumentsProcessed();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.printShackData();

	sweet::Data::Sphere2D::Config sphere2DData_Config;
	sphere2DData_Config.setupAuto(shackSphere2DDataOps);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(sphere2DData_Config);

	LevelSingleton levelSingleton;

	if (shackLibPFASST->nlevels != 1)
	{
		SWEETErrorFatal("For SDC, nlevels has to be equal to 1");
	}

	shackLibPFASST->postprocess_nsweeps();

	// set up LevelSingleton
	levelSingleton.level = 0;
	levelSingleton.sphere2DDataConfig.setupAuto(shackSphere2DDataOps);

	std::cout << "SPH config string: " << levelSingleton.sphere2DDataConfig.getConfigInformationString() << std::endl;

	// setup data operators
	levelSingleton.ops.setup(
		&(levelSingleton.sphere2DDataConfig),
		shackSphere2DDataOps);

	// define the SWEET parameters
	const int nfields = 3; // number of vector fields (here, height and two horizontal velocities)
	int nvars_per_field;
	nvars_per_field = 2 * levelSingleton.sphere2DDataConfig.spectral_array_data_number_of_elements; // number of degrees of freedom per vector field

	// initialize the topography before instantiating the Sphere2DDataCtxSDC object
	/*if (shackBenchmarks->benchmark_name == "flow_over_mountain")
	{
		// create h_topo with the configuration at the finest level
		shackBenchmarks->h_topography = sweet::Data::Sphere2D::DataGrid(&(levelSingleton.sphere2DDataConfig));

		// initialize the topography
		levelSingleton.benchmarks.benchmark->setup_topography();
	}*/

	{
		// instantiate the Sphere2DDataCtx object
		int nnodes[1];
		nnodes[0] = shackLibPFASST->nnodes;
		Sphere2DDataCtxSDC pd_ctx;

		pd_ctx.shackRegistration(&shackProgArgDict);

		pd_ctx.setup(
			&shackProgArgDict,
			&levelSingleton,
			nnodes);

		// get the C string length (needed by Fortran...)
		int string_length = shackLibPFASST->nodes_type.size();

		// flag for the RK stepper
		const int rk_stepper_flag = (shackLibPFASST->use_rk_stepper) ? 1 : 0;

		// call LibPFASST to advance in time
		fmain(
			&pd_ctx,                                             // user defined context
			&shackLibPFASST->niters,                             // number of SDC iterations
			&shackLibPFASST->nsweeps[0],                         // number of SDC sweeps on coarse level
			&nnodes[0],                                          // number of SDC nodes
			(shackLibPFASST->nodes_type).c_str(),                // type of nodes
			&string_length,                                      // length of (shackLibPFASST->nodes_type).c_str()
			&rk_stepper_flag,                                    // flag for the RK stepper => 1 means true, 0 means false
			&nfields,                                            // number of vector fields
			&nvars_per_field,                                    // number of dofs per vector field
			&(pd_ctx.shackTimestepControl->maxSimulationTime), // simulation time
			&(pd_ctx.shackTimestepControl->currentTimestepSize));
	}

	MPI_Finalize();
}
