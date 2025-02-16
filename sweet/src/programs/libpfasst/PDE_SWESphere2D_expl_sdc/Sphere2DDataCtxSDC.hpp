#ifndef PROGRAMS_LIBPFASST_PDE_SWESPHERE2D_EXPL_SDC_SPHERE2DDATACTXSDC_HPP
#define PROGRAMS_LIBPFASST_PDE_SWESPHERE2D_EXPL_SDC_SPHERE2DDATACTXSDC_HPP

#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/IO/Shack.hpp>
#include <vector>
#include "../../PDE_SWESphere2D/Diagnostics.hpp"
#include "../interface/LevelSingleton.hpp"

#include "../../PDE_SWESphere2D/TimeOld/PDESWESphere2DTS_ln_erk.hpp"

#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>

#include "../../PDE_SWESphere2D/Shack.hpp"
#include "../../PDE_SWESphere2D/TimeOld/ShackTimeDiscretization.hpp"
#include "../ShackLibPFASST.hpp"

// Class containing the context necessary to evaluate the right-hand sides
// Currently only contains a pointer to the LevelSingleton and the sweet::Shacks::ShackDictionary object

class Sphere2DDataCtxSDC
{

public:
	sweet::Shacks::Dictionary *shackDict;

	sweet::Data::Sphere2D::Shack *shackSphere2DDataOps;
	sweet::IO::Shack *shackIOData;
	sweet::TimeTree::Shack *shackTimestepControl;
	PDE_SWESphere2D::Shack *shackPDESWESphere2D;
	ShackLibPFASST *shackLibPFASST;
	ShackTimeDiscretization *shackTimeDisc;

	PDE_SWESphere2D::Diagnostics diagnostics;

	bool shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict)
	{
		shackDict = io_shackDict;

		shackSphere2DDataOps = shackDict->getAutoRegistration<sweet::Data::Sphere2D::Shack>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*shackDict);

		shackIOData = shackDict->getAutoRegistration<sweet::IO::Shack>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*shackDict);

		shackPDESWESphere2D = shackDict->getAutoRegistration<PDE_SWESphere2D::Shack>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*shackDict);

		shackTimeDisc = shackDict->getAutoRegistration<ShackTimeDiscretization>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*shackDict);

		shackTimestepControl = shackDict->getAutoRegistration<sweet::TimeTree::Shack>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*shackDict);

		shackLibPFASST = shackDict->getAutoRegistration<ShackLibPFASST>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*shackDict);

		return true;
	}

	Sphere2DDataCtxSDC()
	{
	}

	bool setup(sweet::Shacks::Dictionary *i_shackDict, LevelSingleton *i_singleton, int *i_nnodes)
	{
		shackDict = i_shackDict;
		levelSingleton = i_singleton;

		int rank = 0;
		int nprocs = 0;

		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

		if (!shackDict)
		{
			SWEETErrorFatal("Sphere2DDataCtx: shackDict pointer is NULL!");
		}

		if (!levelSingleton)
		{
			SWEETErrorFatal("Sphere2DDataCtx: levelSingleton pointer is NULL!");
		}

		// initialize the ln_erk time stepper
		timestepper_ln_erk = new PDESWESphere2DTS_ln_erk;
		timestepper_ln_erk->shackRegistration(shackDict);

		// this is never used but this makes clear that with niters=1,
		// we're actually just calling ERK1
		int timestepping_order = 1;
		timestepper_ln_erk->setup(&(*levelSingleton).ops, timestepping_order);

		// initialize the residuals
		residuals.resize(nprocs, std::vector<double>(0, 0.));

		// initialize the diagnostics object
		diagnostics.setup(
			&(*i_singleton).ops,
			shackPDESWESphere2D,
			0);

		return true;
	}

	// Destructor
	~Sphere2DDataCtxSDC()
	{
		delete timestepper_ln_erk;
	}

	// Getter for the sphere2D data configuration
	sweet::Data::Sphere2D::Config *get_sphere2d_data_config() const
	{
		return &(levelSingleton->sphere2DDataConfig);
	}

	// Getter for the sphere2D data configuration
	PDE_SWESphere2D::Benchmarks::BenchmarksCombined *get_swe_benchmark() const
	{
		return &(levelSingleton->benchmarks);
	}

	// Getter for the sphere2D data operators
	sweet::Data::Sphere2D::Operators *get_sphere2d_operators() const
	{
		return &(levelSingleton->ops);
	}

#if 0
	// Getter for the sphere2D data operators with no dealiasing
	sweet::Data::Sphere2D::Operators *get_sphere2d_operators_nodealiasing() const
	{
		return &(levelSingleton->opNoDealiasing);
	}
#endif

	// Getter for the sphere2D diagnostics
	PDE_SWESphere2D::Diagnostics *get_sphere2d_diagnostics()
	{
		return &diagnostics;
	}

	// Getter for the explicit timestepper
	PDESWESphere2DTS_ln_erk *get_ln_erk_timestepper() const
	{
		return timestepper_ln_erk;
	}

	// Getter for the simulationVariables object
	sweet::Shacks::Dictionary *get_simulation_variables() const
	{
		return shackDict;
	}

	// Getter for the number of levels
	int get_number_of_levels() const
	{
		return 1;
	}

	// Save the physical invariants
	void save_grid_invariants(int i_niter)
	{
		time.push_back(shackTimestepControl->currentTimestepSize * i_niter);
		mass.push_back(diagnostics.total_mass);
		energy.push_back(diagnostics.total_energy);
		potentialEnstrophy.push_back(diagnostics.total_potential_enstrophy);
	}

	// Getters for the time and invariants vectors
	const std::vector<double> &get_time() const
	{
		return time;
	}
	const std::vector<double> &get_mass() const { return mass; }
	const std::vector<double> &get_energy() const { return energy; }
	const std::vector<double> &get_potential_enstrophy() const { return potentialEnstrophy; }

	// Getters for the residuals
	const std::vector<std::vector<double>> &get_residuals() const { return residuals; }
	std::vector<std::vector<double>> &get_residuals() { return residuals; }

protected:
	// Pointer to the LevelSingleton object
	LevelSingleton *levelSingleton;

	// Pointer to the ln_erk timestepper
	PDESWESphere2DTS_ln_erk *timestepper_ln_erk;

	// Saved Residuals for each processor
	std::vector<std::vector<double>> residuals;

	Sphere2DDataCtxSDC(const Sphere2DDataCtxSDC &);
	Sphere2DDataCtxSDC &operator=(const Sphere2DDataCtxSDC &);

	// Vectors used for plotting
	std::vector<double> time;
	std::vector<double> mass;
	std::vector<double> energy;
	std::vector<double> potentialEnstrophy;
};

#endif
