/*
 * Author: Francois Hamon & Martin Schreiber <SchreiberX@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/libpfasst/PDE_SWESphere2D_mlsdc/
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
#include "PDE_SWESphere2D_mlsdc/Sphere2DDataCtx.hpp"
#include "../PDE_SWESphere2D/BenchmarksCombined.hpp"
#include "../PDE_SWESphere2D/Diagnostics.hpp"
#include "../PDE_SWESphere2D/TimeOld/ShackTimeDiscretization.hpp"



#define WITH_MPI

extern "C"
{
/* Driver function for pfasst control */
void fmain (
		Sphere2DDataCtx* pd_ctx,
		const int*     nlevels,
		const int*     niters,
		const int      nsweeps[],
		const int      nnodes[],
		const char*    qtype_name,
		const int*     qtype_name_len,
		const int*     use_rk_stepper, // 1 means true, 0 means false
		const int*     nfields,
		const int      nvars_per_field[],
		double*        t_max,
		double*        dt);
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

	PDE_SWESphere2D::Shack *shackPDESWESphere2D = shackProgArgDict.getAutoRegistration<PDE_SWESphere2D::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	ShackLibPFASST *shackLibPFASST = shackProgArgDict.getAutoRegistration<ShackLibPFASST>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	//sweet::TimeTree::ShackTimestepControl *shackTimestepControl =
			shackProgArgDict.getAutoRegistration<sweet::TimeTree::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	//Shack *shackBenchmarks =
			shackProgArgDict.getAutoRegistration<PDE_SWESphere2D::Benchmarks::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	//sweet::IO::ShackIOData *shackIOData =
			shackProgArgDict.getAutoRegistration<sweet::IO::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	//ShackTimeDiscretization *shackTimeDisc =
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


	std::vector<LevelSingleton> levelSingletons;

	// set output time scale to hours
	//shackIOData->output_time_scale = 1.0/(60.0*60.0);
	//shackIOData->output_time_scale_inv = 60.0*60.0;

	if ((shackPDESWESphere2D->viscosity > 0) || (shackPDESWESphere2D->viscosity_order != 2))
	{
		SWEETErrorFatal("To apply viscosity, use the --libpfasst-u2/4/6/8 flags, not -u or -U!");
	}

	shackLibPFASST->postprocess_hyperviscosity();
	shackLibPFASST->postprocess_nsweeps();

	// define the number of levels and SDC nodes for each level
	// note: level #nlevels-1 is the finest, level #0 is the coarsest

	std::vector<int> nnodes;
	nnodes.resize(shackLibPFASST->nlevels);
	nnodes[shackLibPFASST->nlevels-1] = shackLibPFASST->nnodes; // finest level

	switch (shackLibPFASST->nlevels)
	{
	// One level (nothing to do)
	case 1: {
		break;
	}
	// Two levels
	case 2: {
		if (shackLibPFASST->nnodes == 3)
			nnodes[0] = 2;
		else if (shackLibPFASST->nnodes == 5 ||
				shackLibPFASST->nnodes == 9)
			nnodes[0] = 3;
		else if (shackLibPFASST->nnodes == 7)
			nnodes[0] = 4; // for rk_stepper
		else
			SWEETErrorFatal("With 2 levels, the number of SDC nodes on the fine level must be either 3, 5, or 9");
		break;
	}
	// Three levels
	case 3: {
		if (shackLibPFASST->nnodes == 9)
		{
			nnodes[0] = 3;
			nnodes[1] = 5;
		}
		else if (shackLibPFASST->nnodes == 5)
		{
			nnodes[0] = 2;
			nnodes[1] = 3;
		}
		else if (shackLibPFASST->nnodes == 3)
		{
			nnodes[0] = 3;
			nnodes[1] = 3;
		}
		else
			SWEETErrorFatal("With 3 levels, the number of SDC nodes on the fine level must be either 5, or 9");
		break;
	}
	// All other cases not supported yet
	default:
		SWEETErrorFatal("Only 1, 2, or 3 levels are currently supported");
	}

	// setup the LevelSingletons for all levels
	// note: level #nlevels-1 is the finest, level #0 is the coarsest

	levelSingletons.resize(shackLibPFASST->nlevels);

	// setup the finest level singleton
	const int fineLevelId = shackLibPFASST->nlevels-1;


	levelSingletons[fineLevelId].level = fineLevelId;

	// setup data configuration in fine level

	levelSingletons[fineLevelId].sphere2DDataConfig.setupAuto(shackSphere2DDataOps);

	std::cout << "SPH config string: " << levelSingletons[fineLevelId].sphere2DDataConfig.getConfigInformationString() << std::endl;

	// setup data operators in fine level

	levelSingletons[fineLevelId].ops.setup(
			&(levelSingletons[fineLevelId].sphere2DDataConfig),
			shackSphere2DDataOps
	);

	// define the number of modes for the coarser levels
	for (int i = 1; i < shackLibPFASST->nlevels; i++)
	{
		const int thisLevelId = shackLibPFASST->nlevels-1-i;
		levelSingletons[thisLevelId].level = thisLevelId;

		// compute "additional" modes (negative because we're coarsening)
		// use 1 - alpha to compute what to take away (we want to have alpha^i * res modes on level n-1-i)
		double coarsener = 1 - shackLibPFASST->coarsening_multiplier;
		int additional_modes_lat = 1 - std::ceil(shackSphere2DDataOps->space_res_spectral[0]*pow(coarsener,i));
		int additional_modes_lon = 1 - std::ceil(shackSphere2DDataOps->space_res_spectral[1]*pow(coarsener,i));
		// setup data configuration at this level
		levelSingletons[thisLevelId].sphere2DDataConfig.setupAdditionalModes(
				&(levelSingletons[thisLevelId + 1].sphere2DDataConfig),
				additional_modes_lat,
				additional_modes_lon,
				shackSphere2DDataOps
		);

		// setup data operators at this level
		levelSingletons[thisLevelId].ops.setup(
				&(levelSingletons[thisLevelId].sphere2DDataConfig),
				shackSphere2DDataOps
		);
	}

	// define the SWEET parameters

	const int nfields = 3;  // number of vector fields (here, height and two horizontal velocities)
	std::vector<int> nvars_per_field;
	nvars_per_field.resize(shackLibPFASST->nlevels);
	for (int i = 0; i < shackLibPFASST->nlevels; ++i)
	{
		// number of degrees of freedom per vector field
		nvars_per_field[i] = 2*levelSingletons[i].sphere2DDataConfig.spectral_array_data_number_of_elements; 
	}

	{
		// instantiate the Sphere2DDataCtx object
		Sphere2DDataCtx pd_ctx;

		pd_ctx.shackRegistration(&shackProgArgDict);

		pd_ctx.setup(
				&shackProgArgDict,
				&levelSingletons,
				&nnodes
		);

		// get the C string length (needed by Fortran...)
		int string_length = shackLibPFASST->nodes_type.size();

		// flag for the RK stepper
		const int rk_stepper_flag = shackLibPFASST->use_rk_stepper ? 1 : 0;

		// call LibPFASST to advance in time
		fmain(
				&pd_ctx,								// user defined context
				&shackLibPFASST->nlevels,			// number of SDC levels
				&shackLibPFASST->niters,			// number of SDC iterations
				shackLibPFASST->nsweeps.data(),		// number of SDC sweeps on coarse level
				nnodes.data(),						// number of SDC nodes
				(shackLibPFASST->nodes_type).c_str(),	// type of nodes
				&string_length,						// length of (shackLibPFASST->nodes_type).c_str()
				&rk_stepper_flag,					// flag for the RK stepper => 1 means true, 0 means false
				&nfields,							// number of vector fields
				nvars_per_field.data(),				// number of dofs per vector field
				&(pd_ctx.shackTimestepControl->maxSimulationTime),   // simulation time
				&(pd_ctx.shackTimestepControl->currentTimestepSize)  // time step size
		);
	}


	MPI_Finalize();
}

