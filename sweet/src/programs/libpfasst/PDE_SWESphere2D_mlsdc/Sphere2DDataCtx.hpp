#ifndef PROGRAMS_LIBPFASST_PDE_SWESPHERE2D_MLSDC_SPHERE2DDATACTX_HPP
#define PROGRAMS_LIBPFASST_PDE_SWESPHERE2D_MLSDC_SPHERE2DDATACTX_HPP

#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/IO/Shack.hpp>
#include <vector>
#include "../interface/LevelSingleton.hpp"

#include "../../PDE_SWESphere2D/BenchmarksCombined.hpp"
#include "../../PDE_SWESphere2D/Diagnostics.hpp"

#include "../../PDE_SWESphere2D/TimeOld/PDESWESphere2DTS_lg_erk_lc_n_erk.hpp"
#include "../../PDE_SWESphere2D/TimeOld/PDESWESphere2DTS_lg_irk.hpp"


#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>

#include "../../PDE_SWESphere2D/Shack.hpp"
#include "../../PDE_SWESphere2D/TimeOld/ShackTimeDiscretization.hpp"
#include "../ShackLibPFASST.hpp"

/**
 * Class containing the context necessary to evaluate the right-hand sides
 *
 * Currently only contains a pointer to the level singletons and the sweet::Shacks::ShackDictionary object

 */
class Sphere2DDataCtx
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

public:

	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)
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

	Sphere2DDataCtx()
	{
	}

	bool setup(
			sweet::Shacks::Dictionary *i_shackDict,
			std::vector<LevelSingleton> *i_singletons,
			std::vector<int> *i_nnodes
	)
	{
		shackDict = i_shackDict;
		levelSingletons = i_singletons;

		int rank   = 0;
		int nprocs = 0;

		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

		if (!shackDict)
		{
			SWEETErrorFatal("Sphere2DDataCtx: shackDict pointer is NULL!");
		}

		if (!levelSingletons)
		{
			SWEETErrorFatal("Sphere2DDataCtx: levelSingletons pointer is NULL!");
		}

		// use first order integration in time for all pieces (only order supported)
		if (shackTimeDisc->timestepping_order != -1)
		{
			std::cout << "WARNING: Supplying the timestepping order manually is not supported!" << std::endl;
		}

		int timestepping_order = 1;

		// resize vectors
		timestepper_lg_erk_lc_n_erk.resize(levelSingletons->size());
		timestepper_lg_irk.resize(levelSingletons->size());

		// Strang splitting version that should be used:
		int version_id = 1;

		// initialize timesteppers for each level
		for (unsigned int level = 0; level < levelSingletons->size(); level++)
		{
			timestepper_lg_erk_lc_n_erk[level] = new PDESWESphere2DTS_lg_erk_lc_n_erk;
			timestepper_lg_erk_lc_n_erk[level]->shackRegistration(shackDict);
			timestepper_lg_erk_lc_n_erk[level]->setup(&(levelSingletons->at(level).ops), timestepping_order, version_id);

			timestepper_lg_irk[level] = new PDESWESphere2DTS_lg_irk;
			timestepper_lg_irk[level]->shackRegistration(shackDict);
			timestepper_lg_irk[level]->setup(&(levelSingletons->at(level).ops), timestepping_order, shackTimestepControl->currentTimestepSize);
		}

		// initialize the residuals
		residuals.resize(nprocs, std::vector<double>(0,0.));

		diagnostics.setup(
				&((*i_singletons)[shackLibPFASST->nlevels-1].ops),
				shackPDESWESphere2D,
				0
		);

		return true;
	}

	// Destructor
	~Sphere2DDataCtx()
	{
		for (auto & p : timestepper_lg_erk_lc_n_erk)
		{
			delete p;
		}
		for (auto & p : timestepper_lg_irk)
		{
			delete p;
		}
	}

	// Getter for the sphere2D data configuration at level i_level
	sweet::Data::Sphere2D::Config* get_sphere2d_data_config(int i_level) const
	{
		return &(levelSingletons->at(i_level).sphere2DDataConfig);
	}

	// Getter for the sphere2D data configuration at level i_level
	PDE_SWESphere2D::Benchmarks::BenchmarksCombined* get_swe_benchmark(int i_level) const
	{
		return &(levelSingletons->at(i_level).benchmarks);
	}

	// Getter for the sphere2D data operators at level i_level
	sweet::Data::Sphere2D::Operators* get_sphere2d_operators(int i_level) const
	{
		return &(levelSingletons->at(i_level).ops);
	}

#if 0
	// Getter for the sphere2D data operators with no dealiasing at the fine level
	sweet::Data::Sphere2D::Operators* get_sphere2d_operators_nodealiasing() const
	{
		return &(levelSingletons->back().opNoDealiasing);
	}
#endif

	// Getter for the sphere2D diagnostics at the fine level
	PDE_SWESphere2D::Diagnostics* get_sphere2d_diagnostics()
	{
		return &diagnostics;
	}

	// Getter for the explicit timestepper
	PDESWESphere2DTS_lg_erk_lc_n_erk* get_lg_erk_lc_n_erk_timestepper(int i_level) const
	{
		return timestepper_lg_erk_lc_n_erk.at(i_level);
	}

	// Getter for the implicit timestepper
	PDESWESphere2DTS_lg_irk* get_lg_irk_timestepper(int i_level) const
	{
		return timestepper_lg_irk.at(i_level);
	}

	// Getter for the simulationVariables object
	sweet::Shacks::Dictionary* get_simulation_variables() const
	{
		return shackDict;
	}

	// Getter for the number of levels
	int get_number_of_levels() const
	{
		return levelSingletons->size();
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
	const std::vector<double>& get_time()                const { return time; }
	const std::vector<double>& get_mass()                const { return mass; }
	const std::vector<double>& get_energy()              const { return energy; }
	const std::vector<double>& get_potential_enstrophy() const { return potentialEnstrophy; }

	// Getters for the residuals
	const std::vector<std::vector<double> >& get_residuals() const { return residuals; }
	std::vector<std::vector<double> >&       get_residuals()       { return residuals; }

protected:

	// Pointer to the LevelSingletons vector
	std::vector<LevelSingleton> *levelSingletons;

	// Pointer to the lg_erk_lc_n_erk timestepper (used for ceval_f1, ceval_f2)
	std::vector<PDESWESphere2DTS_lg_erk_lc_n_erk*> timestepper_lg_erk_lc_n_erk;

	// Pointer to the lg_irk timestepper (used for ccomp_f2)
	std::vector<PDESWESphere2DTS_lg_irk*> timestepper_lg_irk;

	// Saved Residuals for each processor
	std::vector<std::vector<double> > residuals;


	Sphere2DDataCtx(const Sphere2DDataCtx&);
	Sphere2DDataCtx& operator=(const Sphere2DDataCtx&);

	// Vectors used for plotting
	std::vector<double> time;
	std::vector<double> mass;
	std::vector<double> energy;
	std::vector<double> potentialEnstrophy;

};

#endif
